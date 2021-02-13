#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;
use List::Compare;
use List::Util qw(sum);
use List::Uniq ':all';

my $selfpath=dirname(__FILE__);

my %options;
my @standard_options =("help|h!",
                       "mode=s",
                       "bam=s",
                       "fq=s",
                       "primer_bed=s",
                       "ref=s",
                       "t|thread=s",
                       "o=s",
                       "Rtch=s",
                       "adjleader=s",
                       "adjTRS=s",
                       "covcutoff=s",
                       "extractfasta!",
                       "noguide!"
                       );

GetOptions( \%options, @standard_options );

###################### parameters setting ######################
# if no arguments supplied print the usage and exit
if ((keys (%options))==0) {
    print "please use -h to get help\n";
    exit;
}
# If the -help option is set, print the usage and exit
if ($options{'help'}) {
    print "\nUsage example\:
  perl LeTRS.pl -t 16 -extractfasta -Rtch cDNA -mode noropore -fa example.fastq.gz -primer_bed primer_V3.bed -o LeTRS_output 
  perl LeTRS.pl -t 16 -extractfasta -Rtch RNA -mode noropore -fq example.fastq.gz -primer_bed primer_V3.bed -o LeTRS_output
  perl LeTRS.pl -t 16 -extractfasta -Rtch RNA -mode noropore -fq example.fastq.gz -primer_bed primer_V3.bed -o LeTRS_output -ref reference_folder
  perl LeTRS.pl -t 16 -extractfasta -mode illumia -fq #1.fasq.gz:#2.fasq.gz -primer_bed primer_V3.bed -o LeTRS_output
  perl LeTRS.pl -t 16 -extractfasta -mode illumia -bam example.bam -o LeTRS_output

Required options:
  -mode             \"noropore\" or \"illumia\" of the input fastq file.
  -Rtch             \"RNA\" (direct RNA) or \"cDNA\" (amplicon cDNA) to indicate the sequencing library for \"noropore\" mode.
  -primer_bed       amplicon primer bed file is required for \"noropore cDNA\" and \"illumina\" modes.
  -fq               fastq file (the paired reads can be provide as \"#1.fasq.gz:#2.fasq.gz\").
  -bam              custom bam can also be provided if the \"-fq\" dosen't exist.
  -ref              custom sars-cov-2 or other coronavirus reference folder.

Optional options:
  -extractfasta     to extract the reads contain the identified leader-TRS junctions in fasta format.
  -noguide          this option turns off the known leader-TRS guided alignment.
  -t/-thread        number of threads, 1 by default.
  -o                output path, \"./\" by default.
  -adjleader INT    leader junction boundary tolerance +-10 nts default.
  -adjTRS INT       TRS junction boundary tolerance -20 nts to 1 nt ahead the ATG of knonwn orfs default.
  -covcutoff INT    coverge cutoff for the novel leader-TRS indentification, 10 by default.

  -h/-help          Produce help message.\n\n";
    exit;
}

my $chackbam; my $chackfa; 
# Compulsory items

if (!exists $options{'mode'}) {
    print "mode is requeired\n";
    exit;
}

if (exists ($options{'bam'})) {
    $chackbam=1;
}elsif (!exists ($options{'bam'})){
    $chackbam=0;}
if (exists ($options{'fq'})) {
    $chackfa=1;
}elsif (!exists ($options{'fq'})){
    $chackfa=0;
}
if ($chackbam==0 and $chackfa==0) {
    print "fa or bam is requeired\n";
    exit;
}elsif ($chackbam==1 and $chackfa==0) {
    print "it is running with bam\n";
}elsif ($chackbam==0 and $chackfa==1) {
    print "it is running with fa\n";
}elsif ($chackbam==1 and $chackfa==1) {
    print "please only provide one of bam and fa\n";
}

my $pathtoreference;
if (!exists ($options{'ref'})) {
    $pathtoreference="$selfpath/references_sars_covid_2";
}elsif(exists ($options{'ref'})) {
    $pathtoreference="$options{'ref'}";
}

my $thread;
if (!exists ($options{'t'})) {
    $thread=1;
}elsif (exists ($options{'t'})) {
    $thread= $options{'t'};
}
print "thread: $thread\n";
my $outputpath;
if (!exists ($options{'o'})) {
    $outputpath= "\./";
}elsif (exists ($options{'o'})) {
    $outputpath= $options{'o'};
    mkdir "$outputpath";
}
print "path: $outputpath\n";

###################### parsing ######################
my $leaderadjectnumber;
if (!exists ($options{'adjleader'})) {
    if ($options{'mode'} eq "noropore") {
        $leaderadjectnumber= "10";
    }elsif ($options{'mode'} eq "illumia") {
        $leaderadjectnumber= "10";
    }
}elsif (exists ($options{'adjleader'})) {
    $leaderadjectnumber= $options{'adjleader'};
}
print "junction site at leader adject number: \+\-$leaderadjectnumber\n";

my $orfadjectnumber;
if (!exists ($options{'adjorf'})) {
    if ($options{'mode'} eq "noropore") {
        $orfadjectnumber= "20";
    }elsif ($options{'mode'} eq "illumia") {
        $orfadjectnumber= "20";
    }
}elsif (exists ($options{'adjorf'})) {
    $orfadjectnumber= $options{'adjorf'};
}
print "junction site at orf adject number: \-$orfadjectnumber\n";

my $cutoff=10;
if (!exists ($options{'covcutoff'})) {
    $cutoff= "10";
}elsif (exists ($options{'covcutoff'})) {
    $cutoff= $options{'covcutoff'};
}
print "novel junction cov cutoff: $cutoff\n";

###################### script running ######################
if ($options{'mode'} eq "noropore") {
    print "it is noropore mode now\n";
    if ($chackbam==0 and $chackfa==1) {
        if (!exists ($options{'Rtch'})) {
            print "please indicate the minion seq model: RNA or cDNA in --Rtch\n";
            exit;
        }
        if ($options{'Rtch'} eq "cDNA" && !exists ($options{'primer_bed'})) {
            print "please provide a primer bed file\n";
            exit;
        }
        &mappingnoropore;
        print "looking for the leader-TRS\n";
        &tabfix;
        &bamparse;
        &parseresult;
    }elsif($chackbam==1 and $chackfa==0) {
        if (exists ($options{'primer_bed'})) {
            print "\"-primer_bed\" is not functional with \"-bam\"\n";
            exit;
        }
        $options{'Rtch'}=0;
        &bamnoropore;
        print "looking for the leader-TRS\n";
        &tabfix;
        &bamparse;
        &parseresult;
    }
}elsif ($options{'mode'} eq "illumia") {
    print "it is illumia mode now\n";
    if ($chackbam==0 and $chackfa==1) {
        if (exists ($options{'Rtch'})) {
            print "--Rtch is not functional in illumia mode\n";
            exit;
        }
        $options{'Rtch'}=0;
        if (!exists ($options{'primer_bed'})) {
            print "please provide a primer bed file\n";
            exit;
        }
        &mappingillumia;
        print "looking for the leader-TRS\n";
        &tabfix;
        &bamparse;
        &parseresult;
    }elsif($chackbam==1 and $chackfa==0) {
        if (exists ($options{'primer_bed'})) {
            print "\"-primer_bed\" is not functional with \"-bam\"\n";
            exit;
        }
        $options{'Rtch'}=0;
        &bamillumia;
        print "looking for the leader-TRS\n";
        &tabfix;
        &bamparse;
        &parseresult;
    }
}else {
    print "please select illumia or noropore mode to run\n";
}

###################### Mapping SUBS ######################
sub mappingnoropore {
    my $lable="alignment";
    mkdir "$outputpath/$lable\_output";
    
    if ($options{'Rtch'} eq "cDNA" and exists($options{'noguide'})) {
        print "minion seq model: cDNA\n";
        system ("minimap2 -ax splice -t $thread $pathtoreference/genome.fasta $options{'fq'} | samtools view -@ $thread -q 10 -F 2304 -Sb | samtools sort -@ $thread -o $outputpath/$lable\_output/$lable\.sorted.bam");
    }elsif ($options{'Rtch'} eq "cDNA") {
        print "minion seq model: cDNA\n";
        system ("minimap2 -ax splice --junc-bed $pathtoreference/genome.bed -t $thread $pathtoreference/genome.fasta $options{'fq'} | samtools view -@ $thread -q 10 -F 2304 -Sb | samtools sort -@ $thread -o $outputpath/$lable\_output/$lable\.sorted.bam");
    }elsif ($options{'Rtch'} eq "RNA" and exists($options{'noguide'})) {
        print "minion seq model: RNA\n";
        system ("minimap2 -ax splice -uf -k14 -t $thread $pathtoreference/genome.fasta $options{'fq'} | samtools view -@ $thread -q 10 -F 2320 -Sb | samtools sort -@ $thread -o $outputpath/$lable\_output/$lable\.sorted.bam");
    }elsif ($options{'Rtch'} eq "RNA") {
        print "minion seq model: RNA\n";
        system ("minimap2 -ax splice -uf -k14 --junc-bed $pathtoreference/genome.bed -t $thread $pathtoreference/genome.fasta $options{'fq'} | samtools view -@ $thread -q 10 -F 2320 -Sb | samtools sort -@ $thread -o $outputpath/$lable\_output/$lable\.sorted.bam");
    }elsif ($options{'Rtch'} ne "cDNA" and $options{'Rtch'} ne "RNA") {
        print "please indicate the minion seq model: RNA or cDNA in --Rtch\n";
        exit;
    }

    system ("samtools index $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("portcullis prep -t $thread -o $outputpath/$lable\_output/$lable\_portcullis_out $pathtoreference/genome.fasta  $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("portcullis junc --intron_gff --exon_gff -t $thread --orientation SE -o $outputpath/$lable\_output/$lable\_portcullis_out $outputpath/$lable\_output/$lable\_portcullis_out");
}

sub bamnoropore {
    my $lable="alignment";
    mkdir "$outputpath/$lable\_output";
    system ("samtools view -@ $thread -b -q 10 -F 2304 $options{'bam'} | samtools sort -@ $thread -o $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("portcullis prep -t $thread -o $outputpath/$lable\_output/$lable\_portcullis_out $pathtoreference/genome.fasta $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("portcullis junc --intron_gff --exon_gff -t $thread --orientation SE -o $outputpath/$lable\_output/$lable\_portcullis_out $outputpath/$lable\_output/$lable\_portcullis_out");
}

sub mappingillumia {
    my $lable="alignment";
    mkdir "$outputpath/$lable\_output";
    
    my @fasta=split(/\:/, $options{'fq'});
    if ($#fasta < 1) {
        print "the \"illumia\" mode only supports paired end reads\n";
        exit;
    }
    
    print "$fasta[0]\n";
    print "$fasta[1]\n";
    
    my $indexfileExist = -e "$pathtoreference/genome.1.ht2";
    if ($indexfileExist) {
        print "the hisat2-build has done\n";
    } else {
        system ("hisat2-build -f $pathtoreference/genome.fasta $pathtoreference/genome");
    }
    
    if (exists($options{'noguide'})) {
        system ("hisat2 -p $thread -q -t -x $pathtoreference/genome -1 $fasta[0] -2 $fasta[1] -S $outputpath/$lable\_output/$lable\.sam --summary-file $outputpath/$lable\_output/$lable\.mapping.summary");
    }else{
        system ("hisat2 -p $thread -q -t -x $pathtoreference/genome --known-splicesite-infile $pathtoreference/genome_splicesites.txt -1 $fasta[0] -2 $fasta[1] -S $outputpath/$lable\_output/$lable\.sam --summary-file $outputpath/$lable\_output/$lable\.mapping.summary");
    }
    system ("samtools view -@ $thread -q 10 -F 2316 -Sb $outputpath/$lable\_output/$lable\.sam | samtools sort -@ $thread -o $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("samtools index $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("portcullis prep -t $thread -o $outputpath/$lable\_output/$lable\_portcullis_out $pathtoreference/genome.fasta  $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("portcullis junc --intron_gff --exon_gff -t $thread -o $outputpath/$lable\_output/$lable\_portcullis_out $outputpath/$lable\_output/$lable\_portcullis_out");
}

sub bamillumia {
    my $lable="alignment";
    mkdir "$outputpath/$lable\_output";
    system ("samtools view -@ $thread -b -q 10 -F 2304 $options{'bam'} | samtools sort -@ $thread -o $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("portcullis prep -t $thread -o $outputpath/$lable\_output/$lable\_portcullis_out $pathtoreference/genome.fasta $outputpath/$lable\_output/$lable\.sorted.bam");
    system ("portcullis junc --intron_gff --exon_gff -t $thread -o $outputpath/$lable\_output/$lable\_portcullis_out $outputpath/$lable\_output/$lable\_portcullis_out");
}

sub tabfix {
    open(TAB, "$outputpath/alignment_output/alignment_portcullis_out.junctions.tab");
    open(TABR, ">$outputpath/alignment_output/alignment_portcullis_out.junctions_fixed.tab");
    while (<TAB>) {
        unless (/^\n/){
            if (/^index\trefid/) {
                print TABR;
            }else{
                my @splittab=split(/\t/,$_,10);
                my $endfix=$splittab[5]+2;
                my $leftfix=$splittab[7]+1;
                my $rightfix=$splittab[8]+1;
                print TABR "$splittab[0]\t$splittab[1]\t$splittab[2]\t$splittab[3]\t$splittab[4]\t$endfix\t$splittab[6]\t$leftfix\t$rightfix\t$splittab[9]";
            }
        }
    }
    close TAB; close TABR;
    unlink ("$outputpath/alignment_output/alignment_portcullis_out.junctions.tab");
    rename ("$outputpath/alignment_output/alignment_portcullis_out.junctions_fixed.tab", "$outputpath/alignment_output/alignment_portcullis_out.junctions.tab" );
}


sub bamparse {
    my $inputbam;
    if ($options{'bam'}) {
        $inputbam=$options{'bam'};
    }else{
        $inputbam="$outputpath/alignment_output/alignment.sorted.bam";
    }
    
    system ("samtools view -@ $thread $inputbam | awk '(\$6 ~ /N/)' > $outputpath/alignment_output/alignment.sorted.splice.tmp");
    
    open(BAMEXTRACT, "$outputpath/alignment_output/alignment.sorted.splice.tmp");
    open(BAMEXTRACTR, ">$outputpath/alignment_output/alignment.sorted.splice.tab");
    print BAMEXTRACTR "#ID\tnumberN\tleft\tstart\tend\tright\tjunctionlength\tCIGAR\tread\n";

    while (<BAMEXTRACT>) {
        my @eachcols=split(/\t/);
        my @eachCIGAR=split(/(?<=[A-Z])/,$eachcols[5]);
        my $numberN=grep(/N/,@eachCIGAR);
        my ($idx) = grep { $eachCIGAR[$_] =~/N/ } 0 .. $#eachCIGAR;
        #print $idx;
    
        my $fragment1=0;
        for (my $n=0;$n<$idx;$n++) {
            unless($eachCIGAR[$n]=~/I|S/) {
            $eachCIGAR[$n]=~s/[A-Z]//;
            $fragment1=$fragment1+$eachCIGAR[$n];
            }
        }
        my $start=$eachcols[3]+$fragment1-1;
        $eachCIGAR[$idx]=~s/N//;
        my $end=$start+$eachCIGAR[$idx]+1;
    
        my $fragment2=0;
        for (my $n=$idx+1;$n<$#eachCIGAR+1;$n++) {
            unless($eachCIGAR[$n]=~/I|S/) {
                $eachCIGAR[$n]=~s/[A-Z]//;
                $fragment2=$fragment2+$eachCIGAR[$n];
            }
        }
        my $right=$end+$fragment2-1;

        print BAMEXTRACTR "$eachcols[1]\_\_$eachcols[0]\t$numberN\t$eachcols[3]\t$start\t$end\t$right\t$eachCIGAR[$idx]\t$eachcols[5]\t$eachcols[9]\n";
    }
    close BAMEXTRACT; close BAMEXTRACTR;
    unlink ("$outputpath/alignment_output/alignment.sorted.splice.tmp");
}

my $atgporstionnovel;
my %junctionhash;
sub parseresult {
    mkdir "$outputpath/results";
    open (NOVELDETAILS, ">$outputpath/results/novel_junction_details.tab");
    print NOVELDETAILS "subgenome\tpeak_leader_end\tpeak_TRS_start\tACGAAC\tATG_postion\tknown_ATG\t20_leader_seq\tTRS_seq\tfirst_orf_aa\n";
    open (KNOWNJUNCATIONSDETAILS, ">$outputpath/results/known_junction_details.tab");
    print KNOWNJUNCATIONSDETAILS "subgenome\tpeak_leader_end\tpeak_TRS_start\tACGAAC\tATG_postion\t20_leader_seq\tTRS_seq\tfirst_orf_aa\n";

    open (NOVEL, ">$outputpath/results/novel_junction.tab");
    print NOVEL "subgenome\tleader_end\tTRS_start\tnb_count\tnormalized_count\n";
    
    open (KNOWNJUNCATIONS, ">$outputpath/results/known_junction.tab");
    print KNOWNJUNCATIONS "subgenome\tref_leader_end\tpeak_leader_end\tref_TRS_start\tpeak_TRS_start\tpeak_count\tpeak_normalized_count\tcluster_count\tcluster_normalized_count\n";
    
    open (JBED, "$pathtoreference/junction.bed");
    my @junctions=<JBED>;
    close JBED;
    
    my @mappedcoverage;
    open my $getmappedcoverage, "-|", "samtools flagstat $outputpath/alignment_output/alignment.sorted.bam";
    while (<$getmappedcoverage>) {
        if (/mapped \(/) {
            @mappedcoverage=split(/ /);
            print "number of total mapped read: $mappedcoverage[0]\n";
        }
    }
    close $getmappedcoverage;

    my %hashtab; my %hashtabr; #%hashtab for unkonwn junctions, %hashtabr for known junctions.
    my @allindexcollection;
    open(TAB, "$outputpath/alignment_output/alignment_portcullis_out.junctions.tab");
    while (<TAB>) {
        unless (/^index\trefid/ or /^\n/){
            my @eachtabs=split(/\t/,$_);
            if ($eachtabs[4]<80 && $eachtabs[20]>$cutoff) {
                $hashtab{$eachtabs[0]}="$eachtabs[0]\t$eachtabs[4]\t$eachtabs[5]\t$eachtabs[20]";
                push (@allindexcollection,"$eachtabs[0]");
            }
            if ($eachtabs[4]<266) {
                $hashtabr{$eachtabs[0]}="$eachtabs[0]\t$eachtabs[4]\t$eachtabs[5]\t$eachtabs[20]";
            }
        }
    }
    close TAB;
    
    my @eachleftright=();
    open(BAMEXTRACTR, "$outputpath/alignment_output/alignment.sorted.splice.tab");
    my %hashassp; my %hashasspleftright; my %hashassfa; my %hashpaired; #%hashassp for right site; %hashasspleftright for left and right site; %hashassfa for fasta extaction; %hashpaired for the paired reads
    while (<BAMEXTRACTR>) {
        my @eachassps=split(/\t/);
        if ($eachassps[1] eq 1) {
            my $key; my $value; my $keyleftright; my $valueleftright; my $keyfa; my $valuefa;
            ($key, $value)=("$eachassps[3]\t$eachassps[4]","$eachassps[5]");
            push (@{$hashassp{$key}}, $value);
            
            @eachleftright=($eachassps[0],$eachassps[2],$eachassps[5]);
            ($keyleftright, $valueleftright)=("$eachassps[3]\t$eachassps[4]",[@eachleftright]);
            push (@{$hashasspleftright{$keyleftright}}, $valueleftright);

            ($keyfa, $valuefa)=("$eachassps[3]\t$eachassps[4]","$eachassps[0]\t$eachassps[2]\t$eachassps[5]\t$eachassps[8]");
            push (@{$hashassfa{$keyfa}}, $valuefa);
            
            $hashpaired{$eachassps[0]}="$eachassps[3]\t$eachassps[4]";
        }
    }
    close BAMEXTRACTR;
    
    my @indexcollections;
    my $adjectlederup;
    my $adjectorfup;
    foreach my $junction (@junctions) {
        chomp $junction;
        my @junctionsites=split(/\t/,$junction);
        $junctionhash{$junctionsites[3]}=$junctionsites[0];
        my %hashcluster=();
        
        my @clustercounts=(); #for the normal counts collection 
        my @clusterployAcounts1=(); #for the ployA >= 1 counts collection
        my @clusterployAcounts5=(); #for the ployA >= 5 counts collection
        
        my @clusterleftprimerids=(); #for the left primer ids collection
        my @clusterrightprimerids=(); #for the right primer ids collection
        
        #my @clusterATGsite=(); #for the ATGids collection
        
        #####cluster count of junction
        open(TAB, "$outputpath/alignment_output/alignment_portcullis_out.junctions.tab");
        while (<TAB>) {
            my @eachtabs=split(/\t/);
            unless (/^index\trefid/ or /^\n/){
                $adjectlederup=$junctionsites[1]+$leaderadjectnumber;
                my $adjectlederdown=$junctionsites[1]-$leaderadjectnumber;
                $adjectorfup=$junctionsites[3];
                my $adjectorfdown=$junctionsites[2]-$orfadjectnumber;
                if ($eachtabs[4] < $adjectlederup && $eachtabs[4] > $adjectlederdown && $eachtabs[5] < $adjectorfup && $eachtabs[5] > $adjectorfdown) {
                    push (@indexcollections,$eachtabs[0]);
                    push (@clustercounts,$eachtabs[20]); #cluster uniqe junction count.
                    $hashcluster{$eachtabs[0]}=$eachtabs[20];
                    
                    if ($options{'Rtch'} eq "RNA" && !exists $options{'bam'}) {
                        my $countployAcounts1=grep($_>29870, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        push (@clusterployAcounts1,$countployAcounts1); # to get ployA >= 1 by the sart and end postions in the cluster.
                        my $countployAcounts5=grep($_>29874, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        push (@clusterployAcounts5,$countployAcounts5); # to get ployA >=5 by the sart and end postions in the cluster.
                        #my $clusterATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        #push (@clusterATGsite,$clusterATGsiteid);
                    }
                    
                    if ($options{'Rtch'} eq "cDNA" && !exists $options{'bam'}) {
                        my $countployAcounts1=grep($_>29870, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        push (@clusterployAcounts1,$countployAcounts1); # to get ployA >= 1 by the sart and end postions in the cluster.
                        my $countployAcounts5=grep($_>29874, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        push (@clusterployAcounts5,$countployAcounts5); # to get ployA >=5 by the sart and end postions in the cluster.
                        #my $clusterATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        #push (@clusterATGsite,$clusterATGsiteid);
                        
                        open(PRIMERBED, "$options{'primer_bed'}");
                        while (<PRIMERBED>) {
                            chomp;
                            my @eachprimerbed=split(/\t/);
                            if ($eachprimerbed[2] < $adjectlederup && $eachprimerbed[3]=~/_LEFT/) {
                                my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachtabs[4]\t$eachtabs[5]"}});
                                for (my $n=0; $n<$#leftprimerids+1; $n++) {
                                    push (@clusterleftprimerids,$leftprimerids[$n][0]);
                                }
                            }
                                                        
                            if ($eachprimerbed[1] > 21000 && $eachprimerbed[3]=~/_RIGHT/) {
                                my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachtabs[4]\t$eachtabs[5]"}});
                                for (my $n=0; $n<$#rightprimerids+1; $n++) {
                                    push (@clusterrightprimerids,$rightprimerids[$n][0]);
                                }
                            }                         
                        }
                        close PRIMERBED;
                    }
                    
                    if ($options{'mode'} eq "illumia" && !exists $options{'bam'}) {
                        #my $clusterATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        #push (@clusterATGsite,$clusterATGsiteid);
                        
                        open(PRIMERBED, "$options{'primer_bed'}");
                        while (<PRIMERBED>) {
                            chomp;
                            my @eachprimerbed=split(/\t/);
                            if ($eachprimerbed[2] < $adjectlederup && $eachprimerbed[3]=~/_LEFT/) {   
                                my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachtabs[4]\t$eachtabs[5]"}});
                                for (my $n=0; $n<$#leftprimerids+1; $n++) {
                                    push (@clusterleftprimerids,$leftprimerids[$n][0]);
                                }
                            }
                                                        
                            if ($eachprimerbed[1] > 21000 && $eachprimerbed[3]=~/_RIGHT/) {
                                my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachtabs[4]\t$eachtabs[5]"}});
                                for (my $n=0; $n<$#rightprimerids+1; $n++) {
                                    push (@clusterrightprimerids,$rightprimerids[$n][0]);
                                }
                            }
                        }
                    }
                }
            }
        }
        close TAB;
        
        my $sumclustercounts=sum(@clustercounts);
        my $sumclusterployAcounts1=sum(@clusterployAcounts1);
        my $sumclusterployAcounts5=sum(@clusterployAcounts5);
        #my $sumclusterATGsite=sum(@clusterATGsite);
        my @sumclusterleftprimercounts=uniq(@clusterleftprimerids);
        my $sumclusterleftprimercounts=$#sumclusterleftprimercounts+1;
        my @sumclusterrightprimercounts=uniq(@clusterrightprimerids);
        my $sumclusterrightprimercounts=$#sumclusterrightprimercounts+1;
        
        my $clusterleftrightprimerids = List::Compare->new(\@clusterleftprimerids, \@clusterrightprimerids);
        my @clusterleftrightprimeridsintersection = $clusterleftrightprimerids->get_intersection;
        my @clusterleftrightprimeridsunion = $clusterleftrightprimerids->get_union;
        my $clusterleftrightprimeridsintersectioncount=$#clusterleftrightprimeridsintersection+1;

        #find the paired reads with same junction in "illumia" mode.
        my $numberpairedjucntioncluster=0;
        if ($options{'mode'} eq "illumia") {
            foreach (@clusterleftrightprimeridsunion) {
                my @eachids=split(/\_\_/);
                my $testpair=grep (/\_\_$eachids[1]\b/, @clusterleftrightprimeridsunion);
                if ($testpair == 2) {
                    if ($eachids[0] eq 99) {
                        if ($hashpaired{"99\_\_$eachids[1]"} eq $hashpaired{"147\_\_$eachids[1]"}) {
                           $numberpairedjucntioncluster++;
                        }
                    }
                    if ($eachids[0] eq 163) {
                        if ($hashpaired{"163\_\_$eachids[1]"} eq $hashpaired{"83\_\_$eachids[1]"}) {
                            $numberpairedjucntioncluster++;
                        }
                    }
                }
            }
        }
        $numberpairedjucntioncluster=$numberpairedjucntioncluster*2;
        
        my @clusters= (sort {$hashcluster{$a} <=> $hashcluster{$b}} keys %hashcluster);

        #####peak count of junction
        if ($options{'Rtch'} eq "RNA" && !exists $options{'bam'}) {
            if ($#clusters == -1) {
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t0\t$junctionsites[2]\t0\t0\t0\t0\t0\n";
            }else{
                my @eachtabsagain=split(/\t/,$hashtabr{$clusters[-1]});
                my $covratiopeack=sprintf("%.2f", $eachtabsagain[3]/$mappedcoverage[0]*1000000);
                my $covratiocluster=sprintf("%.2f", $sumclustercounts/$mappedcoverage[0]*1000000);
                #my $covsumclusterATGsite=sprintf("%.2f", $sumclusterATGsite/$mappedcoverage[0]*1000000);
                my $covsumclusterployAcounts1=sprintf("%.2f", $sumclusterployAcounts1/$mappedcoverage[0]*1000000);
                my $covsumclusterployAcounts5=sprintf("%.2f", $sumclusterployAcounts5/$mappedcoverage[0]*1000000);
                
                #my $countATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                my $countployAcounts1=grep($_>29870, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                my $countployAcounts5=grep($_>29874, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                #my $covcountATGsiteid=sprintf("%.2f", $countATGsiteid/$mappedcoverage[0]*1000000);
                my $covcountployAcounts1=sprintf("%.2f", $countployAcounts1/$mappedcoverage[0]*1000000);
                my $covcountployAcounts5=sprintf("%.2f", $countployAcounts5/$mappedcoverage[0]*1000000);

                #my $testlength=$#{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}}+1;
                &extractdetails ($junctionsites[0],$eachtabsagain[1],$eachtabsagain[2]);
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t$eachtabsagain[1]\t$junctionsites[2]\t$eachtabsagain[2]\t$eachtabsagain[3]\($countployAcounts1,$countployAcounts5\)\t$covratiopeack\($covcountployAcounts1,$covcountployAcounts5\)\t$sumclustercounts\($sumclusterployAcounts1,$sumclusterployAcounts5\)\t$covratiocluster\($covsumclusterployAcounts1,$covsumclusterployAcounts5\)\n";
            }
        }elsif ($options{'Rtch'} eq "cDNA" && !exists $options{'bam'}) {
            if ($#clusters == -1) {
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t0\t$junctionsites[2]\t0\t0\t0\t0\t0\n";
            }else{
                my @eachtabsagain=split(/\t/,$hashtabr{$clusters[-1]});
                my $covratiopeack=sprintf("%.2f", $eachtabsagain[3]/$mappedcoverage[0]*1000000);
                my $covratiocluster=sprintf("%.2f", $sumclustercounts/$mappedcoverage[0]*1000000);
                #my $covsumclusterATGsite=sprintf("%.2f", $sumclusterATGsite/$mappedcoverage[0]*1000000);
                my $covsumclusterployAcounts1=sprintf("%.2f", $sumclusterployAcounts1/$mappedcoverage[0]*1000000);
                my $covsumclusterployAcounts5=sprintf("%.2f", $sumclusterployAcounts5/$mappedcoverage[0]*1000000);
                
                #my $countATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                my $countployAcounts1=grep($_>29870, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                my $countployAcounts5=grep($_>29874, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                #my $covcountATGsiteid=sprintf("%.2f", $countATGsiteid/$mappedcoverage[0]*1000000);
                my $covcountployAcounts1=sprintf("%.2f", $countployAcounts1/$mappedcoverage[0]*1000000);
                my $covcountployAcounts5=sprintf("%.2f", $countployAcounts5/$mappedcoverage[0]*1000000);
                
                my @clusterleftprimeridspeak=();
                my @clusterrightprimeridspeak=();
                open(PRIMERBED, "$options{'primer_bed'}");
                while (<PRIMERBED>) {
                    chomp;
                    my @eachprimerbed=split(/\t/);
                    if ($eachprimerbed[2] < $adjectlederup && $eachprimerbed[3]=~/_LEFT/) {
                        my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        for (my $n=0; $n<$#leftprimerids+1; $n++) {
                            push (@clusterleftprimeridspeak,$leftprimerids[$n][0]);
                        }
                    }
           
                    if ($eachprimerbed[1] > 21000 && $eachprimerbed[3]=~/_RIGHT/) {   
                        my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        for (my $n=0; $n<$#rightprimerids+1; $n++) {
                            push (@clusterrightprimeridspeak,$rightprimerids[$n][0]);
                        }
                    }                     
                }
                close PRIMERBED;
                my @sumclusterleftprimercountspeak=uniq(@clusterleftprimeridspeak);
                my $sumclusterleftprimercountspeak=$#sumclusterleftprimercountspeak+1;
                my @sumclusterrightprimercountspeak=uniq(@clusterrightprimeridspeak);
                my $sumclusterrightprimercountspeak=$#sumclusterrightprimercountspeak+1;
                
                
                my $clusterleftrightprimeridspeak = List::Compare->new(\@clusterleftprimeridspeak, \@clusterrightprimeridspeak);
                my @clusterleftrightprimeridsintersectionpeak = $clusterleftrightprimeridspeak->get_intersection;
                my $clusterleftrightprimeridsintersectioncountpeak=$#clusterleftrightprimeridsintersectionpeak+1;
                
                my $covsumclusterleftprimercountspeak=sprintf("%.2f", $sumclusterleftprimercountspeak/$mappedcoverage[0]*1000000);
                my $covsumclusterrightprimercountspeak=sprintf("%.2f", $sumclusterrightprimercountspeak/$mappedcoverage[0]*1000000);
                my $covclusterleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]*1000000);
                
                my $covsumclusterleftprimercounts=sprintf("%.2f", $sumclusterleftprimercounts/$mappedcoverage[0]*1000000);
                my $covsumclusterrightprimercounts=sprintf("%.2f", $sumclusterrightprimercounts/$mappedcoverage[0]*1000000);
                my $covclusterleftrightprimeridsintersectioncount=sprintf("%.2f", $clusterleftrightprimeridsintersectioncount/$mappedcoverage[0]*1000000);
                
                #my $testlength=$#{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}}+1;
                &extractdetails ($junctionsites[0],$eachtabsagain[1],$eachtabsagain[2]);
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t$eachtabsagain[1]\t$junctionsites[2]\t$eachtabsagain[2]\t$eachtabsagain[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$countployAcounts1,$countployAcounts5\)\t$covratiopeack\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covcountployAcounts1,$covcountployAcounts5\)\t$sumclustercounts\($sumclusterleftprimercounts,$sumclusterrightprimercounts,$clusterleftrightprimeridsintersectioncount,$sumclusterployAcounts1,$sumclusterployAcounts5\)\t$covratiocluster\($covsumclusterleftprimercounts,$covsumclusterrightprimercounts,$covclusterleftrightprimeridsintersectioncount,$covsumclusterployAcounts1,$covsumclusterployAcounts5\)\n";
            }
        }elsif ($options{'mode'} eq "illumia" && !exists $options{'bam'}) {
            if ($#clusters == -1) {
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t0\t$junctionsites[2]\t0\t0\t0\t0\t0\n";
            }else{
                my @eachtabsagain=split(/\t/,$hashtabr{$clusters[-1]});
                #my $countATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                #my $covcountATGsiteid=sprintf("%.2f", $countATGsiteid/$mappedcoverage[0]*1000000);
                my $covratiopeack=sprintf("%.2f", $eachtabsagain[3]/$mappedcoverage[0]*1000000);
                my $covratiocluster=sprintf("%.2f", $sumclustercounts/$mappedcoverage[0]*1000000);

                my @clusterleftprimeridspeak=();
                my @clusterrightprimeridspeak=();
                open(PRIMERBED, "$options{'primer_bed'}");
                while (<PRIMERBED>) {
                    chomp;
                    my @eachprimerbed=split(/\t/);
                    if ($eachprimerbed[2] < $adjectlederup && $eachprimerbed[3]=~/_LEFT/) {   
                        my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        for (my $n=0; $n<$#leftprimerids+1; $n++) {
                            push (@clusterleftprimeridspeak,$leftprimerids[$n][0]);
                        }
                    }
                                        
                    if ($eachprimerbed[1] > 21000 && $eachprimerbed[3]=~/_RIGHT/) {
                        my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        for (my $n=0; $n<$#rightprimerids+1; $n++) {
                            push (@clusterrightprimeridspeak,$rightprimerids[$n][0]);
                        }
                    }                   
                }
                close PRIMERBED;
                
                my @sumclusterleftprimercountspeak=uniq(@clusterleftprimeridspeak);
                my $sumclusterleftprimercountspeak=$#sumclusterleftprimercountspeak+1;
                my @sumclusterrightprimercountspeak=uniq(@clusterrightprimeridspeak);
                my $sumclusterrightprimercountspeak=$#sumclusterrightprimercountspeak+1;
         
                
                my $clusterleftrightprimeridspeak = List::Compare->new(\@clusterleftprimeridspeak, \@clusterrightprimeridspeak);
                my @clusterleftrightprimeridsintersectionpeak = $clusterleftrightprimeridspeak->get_intersection;
                my @clusterleftrightprimeridsunionpeak = $clusterleftrightprimeridspeak->get_union;
                my $clusterleftrightprimeridsintersectioncountpeak=$#clusterleftrightprimeridsintersectionpeak+1;
                
                my $covsumclusterleftprimercountspeak=sprintf("%.2f", $sumclusterleftprimercountspeak/$mappedcoverage[0]*1000000);
                my $covsumclusterrightprimercountspeak=sprintf("%.2f", $sumclusterrightprimercountspeak/$mappedcoverage[0]*1000000);
                my $covclusterleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]*1000000);
                
                my $covsumclusterleftprimercounts=sprintf("%.2f", $sumclusterleftprimercounts/$mappedcoverage[0]*1000000);
                my $covsumclusterrightprimercounts=sprintf("%.2f", $sumclusterrightprimercounts/$mappedcoverage[0]*1000000);
                my $covclusterleftrightprimeridsintersectioncount=sprintf("%.2f", $clusterleftrightprimeridsintersectioncount/$mappedcoverage[0]*1000000);
                
                my $numberpairedjucntionpeak=0;
                foreach (@clusterleftrightprimeridsunionpeak) {
                    my @eachids=split(/\_\_/);
                    my $testpair=grep (/\_\_$eachids[1]\b/, @clusterleftrightprimeridsunionpeak);
                    if ($testpair == 2) {
                        if ($eachids[0] eq 99) {
                            if ($hashpaired{"99\_\_$eachids[1]"} eq $hashpaired{"147\_\_$eachids[1]"}) {
                                $numberpairedjucntionpeak++;
                            }
                        }
                        if ($eachids[0] eq 163) {
                            if ($hashpaired{"163\_\_$eachids[1]"} eq $hashpaired{"83\_\_$eachids[1]"}) {
                                $numberpairedjucntionpeak++;
                            }
                        }
                    }
                }
                $numberpairedjucntionpeak=$numberpairedjucntionpeak*2;
                my $covnumberpairedjucntionpeak=sprintf("%.2f", $numberpairedjucntionpeak/$mappedcoverage[0]*1000000);
                my $covnumberpairedjucntioncluster=sprintf("%.2f", $numberpairedjucntioncluster/$mappedcoverage[0]*1000000);
                #my $covsumclusterATGsite=sprintf("%.2f", $sumclusterATGsite/$mappedcoverage[0]*1000000);
                
                #my $testlength=$#{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}}+1;
                &extractdetails ($junctionsites[0],$eachtabsagain[1],$eachtabsagain[2]);
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t$eachtabsagain[1]\t$junctionsites[2]\t$eachtabsagain[2]\t$eachtabsagain[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$numberpairedjucntionpeak\)\t$covratiopeack\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covnumberpairedjucntionpeak\)\t$sumclustercounts\($sumclusterleftprimercounts,$sumclusterrightprimercounts,$clusterleftrightprimeridsintersectioncount,$numberpairedjucntioncluster\)\t$covratiocluster\($covsumclusterleftprimercounts,$covsumclusterrightprimercounts,$covclusterleftrightprimeridsintersectioncount,$covnumberpairedjucntioncluster\)\n";
            }
        }else {
            if ($#clusters == -1) {
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t0\t$junctionsites[2]\t0\t0\t0\t0\t0\n";
            }else{
                my @eachtabsagain=split(/\t/,$hashtabr{$clusters[-1]});
                my $covratiopeack=sprintf("%.2f", $eachtabsagain[3]/$mappedcoverage[0]*1000000);
                my $covratiocluster=sprintf("%.2f", $sumclustercounts/$mappedcoverage[0]*1000000);
                
                &extractdetails ($junctionsites[0],$eachtabsagain[1],$eachtabsagain[2]);
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t$eachtabsagain[1]\t$junctionsites[2]\t$eachtabsagain[2]\t$eachtabsagain[3]\t$covratiopeack\t$sumclustercounts\t$covratiocluster\n";
            }
        }
        
        ######extarct known fasta from bam
        if ($options{'extractfasta'}) {
            mkdir "$outputpath/results/fasta";
            open(GETFATSTA, ">$outputpath/results/fasta/leader-$junctionsites[0]\.fa");
            unless ($#clusters == -1) {
                my @eachtabsagain=split(/\t/,$hashtabr{$clusters[-1]});
                foreach (@{$hashassfa{"$eachtabsagain[1]\t$eachtabsagain[2]"}}) {
                    my @getfastas=split(/\t/);
                    print GETFATSTA "\>$getfastas[0] left\:$getfastas[1] start\:$eachtabsagain[1] end\:$eachtabsagain[2] right\:$getfastas[2]\n";
                    print GETFATSTA "$getfastas[3]\n";
                }
            }
            close GETFATSTA;
        }
    }


    my $indexnovel = List::Compare->new(\@indexcollections, \@allindexcollection);
    
    my @indexnovelcollections = $indexnovel->get_complement;
    
    if ($options{'Rtch'} eq "RNA" && !exists $options{'bam'}) {
        my $n=0;
        foreach my $indexnovelcollection(@indexnovelcollections) {
            $n++;
            my @eachnvelintable=split(/\t/,"$hashtab{$indexnovelcollection}");
            my $countployAcounts1novel=grep($_>29870, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            my $countployAcounts5novel=grep($_>29874, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            #my $countATGsiteid=grep($_>$atgporstionnovel+1, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            my $covnoveluniqe=sprintf("%.2f", $eachnvelintable[3]/$mappedcoverage[0]*1000000);
            my $covcountployAcounts1novel=sprintf("%.2f", $countployAcounts1novel/$mappedcoverage[0]*1000000);
            my $covcountployAcounts5novel=sprintf("%.2f", $countployAcounts5novel/$mappedcoverage[0]*1000000);
            #my $covcountATGsiteid=sprintf("%.2f", $countATGsiteid/$mappedcoverage[0]*1000000);
            #my $testlength=$#{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}}+1;
            
            &extractdetailsnovel($n,$eachnvelintable[1],$eachnvelintable[2]);
            print NOVEL "$n\t$eachnvelintable[1]\t$eachnvelintable[2]\t$eachnvelintable[3]\($countployAcounts1novel,$countployAcounts5novel\)\t$covnoveluniqe\($covcountployAcounts1novel,$covcountployAcounts5novel\)\n";
        }
    }elsif ($options{'Rtch'} eq "cDNA" && !exists $options{'bam'}) {
        my $n=0;
        foreach my $indexnovelcollection(@indexnovelcollections) {
            $n++;
            my @eachnvelintable=split(/\t/,"$hashtab{$indexnovelcollection}");
            my $countployAcounts1novel=grep($_>29870, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            my $countployAcounts5novel=grep($_>29874, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            #my $countATGsiteid=grep($_>$atgporstionnovel+1, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            my $covnoveluniqe=sprintf("%.2f", $eachnvelintable[3]/$mappedcoverage[0]*1000000);
            my $covcountployAcounts1novel=sprintf("%.2f", $countployAcounts1novel/$mappedcoverage[0]*1000000);
            my $covcountployAcounts5novel=sprintf("%.2f", $countployAcounts5novel/$mappedcoverage[0]*1000000);
            #my $covcountATGsiteid=sprintf("%.2f", $countATGsiteid/$mappedcoverage[0]*1000000);
            
            my @clusterleftprimeridspeak=();
            my @clusterrightprimeridspeak=();
            open(PRIMERBED, "$options{'primer_bed'}");
            while (<PRIMERBED>) {
                chomp;
                my @eachprimerbed=split(/\t/);
                if ($eachprimerbed[2] < $adjectlederup && $eachprimerbed[3]=~/_LEFT/) {    
                    my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    for (my $n=0; $n<$#leftprimerids+1; $n++) {
                        push (@clusterleftprimeridspeak,$leftprimerids[$n][0]);
                    }
                }
                                
                if ($eachprimerbed[1] > $adjectlederup && $eachprimerbed[3]=~/_RIGHT/) {     
                    my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    for (my $n=0; $n<$#rightprimerids+1; $n++) {
                        push (@clusterrightprimeridspeak,$rightprimerids[$n][0]);
                    }
                }                
            }
            close PRIMERBED;

            my @sumclusterleftprimercountspeak=uniq(@clusterleftprimeridspeak);
            my $sumclusterleftprimercountspeak=$#sumclusterleftprimercountspeak+1;
            my @sumclusterrightprimercountspeak=uniq(@clusterrightprimeridspeak);
            my $sumclusterrightprimercountspeak=$#sumclusterrightprimercountspeak+1;
                
            my $clusterleftrightprimeridspeak = List::Compare->new(\@clusterleftprimeridspeak, \@clusterrightprimeridspeak);
            my @clusterleftrightprimeridsintersectionpeak = $clusterleftrightprimeridspeak->get_intersection;
            my $clusterleftrightprimeridsintersectioncountpeak=$#clusterleftrightprimeridsintersectionpeak+1;
                
            my $covsumclusterleftprimercountspeak=sprintf("%.2f", $sumclusterleftprimercountspeak/$mappedcoverage[0]*1000000);
            my $covsumclusterrightprimercountspeak=sprintf("%.2f", $sumclusterrightprimercountspeak/$mappedcoverage[0]*1000000);
            my $covclusterleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]*1000000);
            
            #my $testlength=$#{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}}+1;
            
            &extractdetailsnovel($n,$eachnvelintable[1],$eachnvelintable[2]);
            print NOVEL "$n\t$eachnvelintable[1]\t$eachnvelintable[2]\t$eachnvelintable[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$countployAcounts1novel,$countployAcounts5novel\)\t$covnoveluniqe\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covcountployAcounts1novel,$covcountployAcounts5novel\)\n";
        }
    }elsif ($options{'mode'} eq "illumia" && !exists $options{'bam'}) {
        my $n=0;
        foreach my $indexnovelcollection(@indexnovelcollections) {
            $n++;
            my @eachnvelintable=split(/\t/,"$hashtab{$indexnovelcollection}");
            #my $countATGsiteid=grep($_>$atgporstionnovel+1, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            #my $covcountATGsiteid=sprintf("%.2f", $countATGsiteid/$mappedcoverage[0]*1000000);
            my $covnoveluniqe=sprintf("%.2f", $eachnvelintable[3]/$mappedcoverage[0]*1000000);

            my @clusterleftprimeridspeak=();
            my @clusterrightprimeridspeak=();
            open(PRIMERBED, "$options{'primer_bed'}");
            while (<PRIMERBED>) {
                chomp;
                my @eachprimerbed=split(/\t/);
                if ($eachprimerbed[2] < $adjectlederup && $eachprimerbed[3]=~/_LEFT/) {       
                    my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    for (my $n=0; $n<$#leftprimerids+1; $n++) {
                          push (@clusterleftprimeridspeak,$leftprimerids[$n][0]);
                    }
                }
                                
                if ($eachprimerbed[1] > $adjectlederup && $eachprimerbed[3]=~/_RIGHT/) {     
                    my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    for (my $n=0; $n<$#rightprimerids+1; $n++) {
                        push (@clusterrightprimeridspeak,$rightprimerids[$n][0]);
                    }
                }                           
            }
            close PRIMERBED;

            my @sumclusterleftprimercountspeak=uniq(@clusterleftprimeridspeak);
            my $sumclusterleftprimercountspeak=$#sumclusterleftprimercountspeak+1;
            my @sumclusterrightprimercountspeak=uniq(@clusterrightprimeridspeak);
            my $sumclusterrightprimercountspeak=$#sumclusterrightprimercountspeak+1;
                
            my $clusterleftrightprimeridspeak = List::Compare->new(\@clusterleftprimeridspeak, \@clusterrightprimeridspeak);
            my @clusterleftrightprimeridsintersectionpeak = $clusterleftrightprimeridspeak->get_intersection;
            my @clusterleftrightprimeridsunionpeak = $clusterleftrightprimeridspeak->get_union;
            my $clusterleftrightprimeridsintersectioncountpeak=$#clusterleftrightprimeridsintersectionpeak+1;
                
            my $covsumclusterleftprimercountspeak=sprintf("%.2f", $sumclusterleftprimercountspeak/$mappedcoverage[0]*1000000);
            my $covsumclusterrightprimercountspeak=sprintf("%.2f", $sumclusterrightprimercountspeak/$mappedcoverage[0]*1000000);
            my $covclusterleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]*1000000);
            
            my $numberpairedjucntionpeak=0;
            foreach (@clusterleftrightprimeridsunionpeak) {
                my @eachids=split(/\_\_/);
                my $testpair=grep (/\_\_$eachids[1]\b/, @clusterleftrightprimeridsunionpeak);
                if ($testpair == 2) {
                    if ($eachids[0] eq 99) {
                        if ($hashpaired{"99\_\_$eachids[1]"} eq $hashpaired{"147\_\_$eachids[1]"}) {
                        $numberpairedjucntionpeak++;
                        }
                    }
                    if ($eachids[0] eq 163) {
                        if ($hashpaired{"163\_\_$eachids[1]"} eq $hashpaired{"83\_\_$eachids[1]"}) {
                            $numberpairedjucntionpeak++;
                        }
                    }
                }
            }
            $numberpairedjucntionpeak=$numberpairedjucntionpeak*2;
            my $covnumberpairedjucntionpeak=sprintf("%.2f", $numberpairedjucntionpeak/$mappedcoverage[0]*1000000);
                
            #my $testlength=$#{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}}+1;
            &extractdetailsnovel($n,$eachnvelintable[1],$eachnvelintable[2]);
            print NOVEL "$n\t$eachnvelintable[1]\t$eachnvelintable[2]\t$eachnvelintable[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$numberpairedjucntionpeak\)\t$covnoveluniqe\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covnumberpairedjucntionpeak\)\n";
        }
    }else{
        my $n=0;
        foreach my $indexnovelcollection(@indexnovelcollections) {
            $n++;
            my @eachnvelintable=split(/\t/,"$hashtab{$indexnovelcollection}");
            #my $countATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            my $covnoveluniqe=sprintf("%.2f", $eachnvelintable[3]/$mappedcoverage[0]*1000000);
            #my $testlength=$#{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}}+1;
            &extractdetailsnovel($n,$eachnvelintable[1],$eachnvelintable[2]);
            print NOVEL "$n\t$eachnvelintable[1]\t$eachnvelintable[2]\t$eachnvelintable[3]\t$covnoveluniqe\n";
        }        
    }
    
    ######extarct unknown fasta from bam
    if ($options{'extractfasta'}) {
        print "extracting fasta with leader-TRS\n";
        mkdir "$outputpath/results/fasta_novel";
        my $n=0;
        foreach my $indexnovelcollection(@indexnovelcollections) {
            $n++;
            #print "$n\n";
            open(GETFATSTANOVEL, ">$outputpath/results/fasta_novel/leader-TRS_novel_$n\.fa");
            my @eachnvelintable=split(/\t/,"$hashtab{$indexnovelcollection}");
            foreach (@{$hashassfa{"$eachnvelintable[1]\t$eachnvelintable[2]"}}) {
                my @getfastas=split(/\t/);
                print GETFATSTANOVEL "\>$getfastas[0] left\:$getfastas[1] start\:$eachnvelintable[1] end\:$eachnvelintable[2] right\:$getfastas[2]\n";
                print GETFATSTANOVEL "$getfastas[3]\n";
            }
            close GETFATSTANOVEL;
        }                    
    }
    
    ######export result
    if ($options{'Rtch'} eq "RNA" && !exists $options{'bam'}) {
        print KNOWNJUNCATIONS "The numbers in the bracket are (reads with > 1 poly A, reads with > 5 poly A)\.\n";
        print KNOWNJUNCATIONS "Normalized count=(Read count/Total number of read mapped on reference genome)*1000000\.\n";
        print KNOWNJUNCATIONS "Total number of read mapped on reference genome is $mappedcoverage[0], excluding mapped reads on reverse strand, not primary alignment and supplementary alignment\.\n";
        print NOVEL "The numbers in the bracket are (reads with > 1 poly A, reads with > 5 poly A)\.\n";
        print NOVEL "Normalized count=(Read count/Total number of read mapped on reference genome)*1000000\.\n";
        print NOVEL "Total number of read mapped on reference genome is $mappedcoverage[0], excluding mapped reads on reverse strand, not primary alignment and supplementary alignment\.\n";
    }elsif ($options{'Rtch'} eq "cDNA" && !exists $options{'bam'}) {
        print KNOWNJUNCATIONS "The numbers in the bracket are (reads with left primers, reads with right primers, reads with both primers, reads with > 1 poly A, reads with > 5 poly A)\.\n";
        print KNOWNJUNCATIONS "Normalized count=(Read count/Total number of read mapped on reference genome)*1000000\.\n";
        print KNOWNJUNCATIONS "Total number of read mapped on reference genome is $mappedcoverage[0], excluding the mapped reads not primary alignment and supplementary alignment\.\n";
        print NOVEL "The numbers in the bracket are (reads with left primers, reads with right primers, reads with both primers, reads with > 1 poly A, reads with > 5 poly A)\.\n";
        print NOVEL "Normalized count=(Read count/Total number of read mapped on reference genome)*1000000\.\n";
        print NOVEL "Total number of read mapped on reference genome doesn't include the mapped reads not primary alignment and supplementary alignment\.\n";
    }elsif ($options{'mode'} eq "illumia" && !exists $options{'bam'}) {
        print KNOWNJUNCATIONS "The numbers in the bracket are (reads with left primers, reads with right primers, reads with both primers, same junction on paired reads with at least a primer)\.\n";
        print KNOWNJUNCATIONS "Normalized count=(Read count/Total number of read mapped on reference genome)*1000000\.\n";
        print KNOWNJUNCATIONS "Total number of read mapped on reference genome is $mappedcoverage[0], excluding the mapped reads unpaired, not primary alignment and supplementary alignment\.\n";
        print NOVEL "The numbers in the bracket are (reads with left primers, reads with right primers, reads with both primers, same junction on paired reads with at least a primer)\.\n";
        print NOVEL "Normalized count=(Read count/Total number of read mapped on reference genome)*1000000\.\n";
        print NOVEL "Total number of read mapped on reference is $mappedcoverage[0], excluding the mapped reads unpaired, not primary alignment and supplementary alignment\.\n";
    }
    close NOVEL; close KNOWNJUNCATIONS;
    close KNOWNJUNCATIONSDETAILS; close NOVELDETAILS;
}

print "please find the leader-TRSs at results folder\n\n";

######extarct details
sub extractdetails {
    my ($each0,$each2,$each4)=($_[0],$_[1],$_[2]);
    #See Bio::Tools::CodonTable for information on the format of a codon table.

    my $inputfilename="$pathtoreference/genome.fasta";
    if ($each2==0) {
        #code
    }else{
    print KNOWNJUNCATIONSDETAILS "$each0\t$each2\t$each4\t";
    my $in = Bio::SeqIO->new(-file => "$inputfilename", -format => 'Fasta');
        while ( my $seq = $in->next_seq() ) {
            my $genomelength=$seq->length;
            my $leaderseq=$seq->subseq($each2-19,$each2);
            my $TRSseqall=$seq->subseq($each4,$genomelength);
            $TRSseqall =~ m/ATG/ig;
            my $atgporstion=pos($TRSseqall)+$each4-3;
            my $TRSseq;
            if ($each4 < $atgporstion) {
                $TRSseq=$seq->subseq($each4,$atgporstion-1);
            }else {
                $TRSseq="\-";
            }
    
            if ($TRSseq =~ m/ACGAAC/ig) {
                print KNOWNJUNCATIONSDETAILS "yes\t";
            }else {
                print KNOWNJUNCATIONSDETAILS "no\t";
            }

            my $prot_obj =$seq->trunc($each4,$genomelength)->translate(-orf => 1, -start => "atg");
            print KNOWNJUNCATIONSDETAILS "$atgporstion\t";
            print KNOWNJUNCATIONSDETAILS "$leaderseq\t";
            print KNOWNJUNCATIONSDETAILS "$TRSseq\t";
            print KNOWNJUNCATIONSDETAILS $prot_obj->seq, "\n";
        }
    }
}
    
sub extractdetailsnovel {
    my $inputfilename="$pathtoreference/genome.fasta";
    my ($each0,$each1,$each2)=($_[0],$_[1],$_[2]);
    print NOVELDETAILS "$each0\t$each1\t$each2\t";
            
    my $in = Bio::SeqIO->new(-file => "$inputfilename", -format => 'Fasta');
    while ( my $seq = $in->next_seq() ) {
        my $genomelength=$seq->length;
        my $leaderseq=$seq->subseq($each1-19,$each1);
        my $TRSseqall=$seq->subseq($each2,$genomelength);            
        $TRSseqall =~ m/ATG/ig;
        my $atgporstion=pos($TRSseqall)+$each2-3;
        $atgporstionnovel=$atgporstion;
        
        my $known_ATG;
        if (exists $junctionhash{$atgporstion}) {
            $known_ATG=$junctionhash{$atgporstion};
        }else{
            $known_ATG="\-";
        }
        
        my $TRSseq;
        if ($each2 < $atgporstion) {
            $TRSseq=$seq->subseq($each2,$atgporstion-1);
        }else {
            $TRSseq="\-";
        }

        if ($TRSseq =~ m/ACGAAC/ig) {
            print NOVELDETAILS "yes\t";
        }else {
            print NOVELDETAILS "no\t";
        }

        
        my $prot_obj =$seq->trunc($each2,$genomelength)->translate(-orf => 1, -start => "atg");
        print NOVELDETAILS "$atgporstion\t";
        print NOVELDETAILS "$known_ATG\t";
        print NOVELDETAILS "$leaderseq\t";
        print NOVELDETAILS "$TRSseq\t";
        print NOVELDETAILS $prot_obj->seq, "\n";
    }
}
