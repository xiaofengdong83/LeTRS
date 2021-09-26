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
                       "version|v!",
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
                       "poscutoff=s",
                       "covcutoff=s",
                       "extractfasta!",
                       "noguide!",
                       "ployALoc=s",
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
  perl LeTRS.pl -t 16 -extractfasta -Rtch cDNA -mode nanopore -fq example.fastq.gz -primer_bed path_to_primer_V3.bed -o LeTRS_output 
  perl LeTRS.pl -t 16 -extractfasta -Rtch RNA -mode nanopore -fq example.fastq.gz -o LeTRS_output
  perl LeTRS.pl -t 16 -extractfasta -Rtch cDNA -mode nanopore -fq example.fastq.gz -primer_bed path_to_custom_primer.bed -o LeTRS_output -ref reference_folder
  perl LeTRS.pl -t 16 -extractfasta -mode illumia -fq #1.fastq.gz:#2.fastq.gz -primer_bed path_to_primer_V3.bed -o LeTRS_output
  perl LeTRS.pl -t 16 -extractfasta -mode illumia -bam example.bam -o LeTRS_output

Required options:
  -mode             \"nanopore\" or \"illumia\" of the input fastq file.
  -Rtch             \"RNA\" (direct RNA) or \"cDNA\" (amplicon cDNA) to indicate the sequencing library for \"nanopore\" mode.
  -primer_bed       amplicon primer bed file is required for \"nanopore cDNA\" and \"illumina\" modes.
  -fq               fastq file (the paired reads can be provide as \"#1.fastq.gz:#2.fastq.gz\").

Optional options:
  -bam              custom bam can also be provided if the \"-fq\" dosen't exist.
  -ref              custom sars-cov-2 or other coronavirus reference folder.
  -extractfasta     to extract the reads contain the identified leader-TRS junctions in fasta format.
  -noguide          this option turns off the known leader-TRS guided alignment.
  -t/-thread        number of threads, 1 by default.
  -o                output path, \"./\" by default.
  -adjleader INT    leader junction boundary tolerance, +-10 nts by default.
  -adjTRS INT       TRS junction boundary tolerance, -20 nts to 1 nt ahead the ATG of knonwn orfs by default.
  -poscutoff INT    postion of leader end cutoff for the novel leader-TRS identification, 80 by default.
  -covcutoff INT    coverge cutoff for the novel leader-TRS identification, 10 by default.
  -ployALoc INT     The position of the first A in ployA tail on the virus genome, auto detection by default.
  
  -v/-version       Print version information.
  -h/-help          Print help message.\n\n";
    exit;
}

if ($options{'version'}) {
    print "v1.0.1\n";
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

my $polyapostion;
my $polyapostion5;
if (!exists ($options{'ployALoc'})) {
    my $inputfilename="$pathtoreference/genome.fasta";
    my $in = Bio::SeqIO->new(-file => "$inputfilename", -format => 'Fasta');
    while ( my $seq = $in->next_seq() ) {
        my $wholegenomeseq= $seq->seq;
        my $wholegenomelength= $seq->length;
        my @eachwholegenomeseq=split(//,$wholegenomeseq);
        my $stepseq=0;
        for (;;) {
            $stepseq--;
            last if ($eachwholegenomeseq[$stepseq] ne "A");
        }
        $wholegenomelength=$wholegenomelength+$stepseq+2;
        if ($stepseq == -1) {
            print "please add ployA sequences in your reference genome";
            exit;
        }else {
            $polyapostion=$wholegenomelength;
            $polyapostion5=$polyapostion+4;
            print "the position of the first A in ployA tail is $polyapostion on the virus genome.\n";
        }
    }
}elsif (exists ($options{'ployALoc'})) {
    $polyapostion= $options{'ployALoc'};
    $polyapostion5=$polyapostion+4;
    print "the position of the first A in ployA tail is $polyapostion on the virus genome.\n";
}


###################### parsing ######################
my $leaderadjectnumber;
if (!exists ($options{'adjleader'})) {
    if ($options{'mode'} eq "nanopore") {
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
    if ($options{'mode'} eq "nanopore") {
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
print "novel junction covervage cutoff: $cutoff\n";

my $cutoffpos=80;
if (!exists ($options{'poscutoff'})) {
    $cutoffpos= "80";
}elsif (exists ($options{'poscutoff'})) {
    $cutoffpos= $options{'poscutoff'};
}
print "novel junction of genome position of leader end cutoff: $cutoff\n";

###################### script running ######################
if ($options{'mode'} eq "nanopore") {
    print "it is nanopore mode now\n";
    if ($chackbam==0 and $chackfa==1) {
        if (!exists ($options{'Rtch'})) {
            print "please indicate the minion seq model: RNA or cDNA in --Rtch\n";
            exit;
        }
        if ($options{'Rtch'} eq "cDNA" && !exists ($options{'primer_bed'})) {
            print "please provide a primer bed file\n";
            exit;
        }
        &mappingnanopore;
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
        &bamnanopore;
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
        &bamparseillumina;
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
    print "please select illumia or nanopore mode to run\n";
}

###################### Mapping SUBS ######################
sub mappingnanopore {
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

sub bamnanopore {
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
    
    my @bamextracts=<BAMEXTRACT>;
    foreach (@bamextracts) {
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

sub bamparseillumina {
    my $inputbam;
    if ($options{'bam'}) {
        $inputbam=$options{'bam'};
    }else{
        $inputbam="$outputpath/alignment_output/alignment.sorted.bam";
    }
    
    system ("samtools view -@ $thread $inputbam > $outputpath/alignment_output/alignment.sorted.splice.tmp");
    
    open(BAMEXTRACT, "$outputpath/alignment_output/alignment.sorted.splice.tmp");
    open(BAMEXTRACTR, ">$outputpath/alignment_output/alignment.sorted.splice.tab");
    print BAMEXTRACTR "#ID1\tnumberN1\tleft1\tstart1\tend1\tright1\tjunctionlength1\tCIGAR1\tread1\t#ID2\tnumberN2\tleft2\tstart2\tend2\tright2\tjunctionlength2\tCIGAR2\tread2\n";

    my %hashbamextracts=();
    my @bamextracts=<BAMEXTRACT>;
    foreach (@bamextracts) {
        chomp;
        my @eachbamextracts=split(/\t/);
        my @eachCIGAR=split(/(?<=[A-Z])/,$eachbamextracts[5]);
        my $numberN=grep(/N/,@eachCIGAR);
        if ($numberN == 0 or $numberN == 1) {
            my @groupbamextracts=("$eachbamextracts[0]","$eachbamextracts[1]","$eachbamextracts[3]","$eachbamextracts[5]","$eachbamextracts[9]");
            push (@{$hashbamextracts{$eachbamextracts[0]}}, [@groupbamextracts]);
        }
    }
    
    for my $keybamextracts (keys %hashbamextracts) {
        #print "$keybamextracts => ";
        if ($#{$hashbamextracts{$keybamextracts}} == 1) {
            my @numgroups;
            if ($hashbamextracts{$keybamextracts}[0][2] <= $hashbamextracts{$keybamextracts}[1][2]) {
                @numgroups=("0","1");
            }
            if ($hashbamextracts{$keybamextracts}[0][2] > $hashbamextracts{$keybamextracts}[1][2]) {
                @numgroups=("1","0");
            }
            
            my @eachCIGARgroup1=split(/(?<=[A-Z])/,$hashbamextracts{$keybamextracts}[0][3]);
            my @eachCIGARgroup2=split(/(?<=[A-Z])/,$hashbamextracts{$keybamextracts}[1][3]);
            my $numberNgroup1=grep(/N/,@eachCIGARgroup1);
            my $numberNgroup2=grep(/N/,@eachCIGARgroup2);
            
            
            if ($numberNgroup1 == 1 or $numberNgroup2 ==1) {
                foreach my $numgroup(@numgroups) {
                    #print "$hashbamextracts{$keybamextracts}[$numgroup][1]\__$hashbamextracts{$keybamextracts}[$numgroup][0]", "\n";
                    my @eachCIGAR=split(/(?<=[A-Z])/,$hashbamextracts{$keybamextracts}[$numgroup][3]);
                    my $numberN=grep(/N/,@eachCIGAR);
        
                    if ($numberN == 1) {
                        my ($idx) = grep { $eachCIGAR[$_] =~/N/ } 0 .. $#eachCIGAR;
                        #print $idx;
    
                        my $fragment1=0;
                        for (my $n=0;$n<$idx;$n++) {
                            unless($eachCIGAR[$n]=~/I|S/) {
                            $eachCIGAR[$n]=~s/[A-Z]//;
                            $fragment1=$fragment1+$eachCIGAR[$n];
                            }
                        }
                        my $start=$hashbamextracts{$keybamextracts}[$numgroup][2]+$fragment1-1;
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
                        print BAMEXTRACTR "$hashbamextracts{$keybamextracts}[$numgroup][1]\_\_$hashbamextracts{$keybamextracts}[$numgroup][0]\t$numberN\t$hashbamextracts{$keybamextracts}[$numgroup][2]\t$start\t$end\t$right\t$eachCIGAR[$idx]\t$hashbamextracts{$keybamextracts}[$numgroup][3]\t$hashbamextracts{$keybamextracts}[$numgroup][4]\t";
                    }
                    if ($numberN == 0) {
                        my $fragment0=0;
                        for (my $n=0;$n<$#eachCIGAR+1;$n++) {
                            unless($eachCIGAR[$n]=~/I|S/) {
                                $eachCIGAR[$n]=~s/[A-Z]//;
                                $fragment0=$fragment0+$eachCIGAR[$n];
                            }
                        }
                        my $right=$hashbamextracts{$keybamextracts}[$numgroup][2]+$fragment0-1;
                        print BAMEXTRACTR "$hashbamextracts{$keybamextracts}[$numgroup][1]\_\_$hashbamextracts{$keybamextracts}[$numgroup][0]\t$numberN\t$hashbamextracts{$keybamextracts}[$numgroup][2]\t\-\t\-\t$right\t\-\t$hashbamextracts{$keybamextracts}[$numgroup][3]\t$hashbamextracts{$keybamextracts}[$numgroup][4]\t";
                    }
                }
                print BAMEXTRACTR "\n";
            }        
        }
    }
    close BAMEXTRACT; close BAMEXTRACTR;
    unlink ("$outputpath/alignment_output/alignment.sorted.splice.tmp");
}

my $atgporstionnovel;
my %junctionhash;
sub parseresult {
    mkdir "$outputpath/results";
    open (NOVELDETAILS, ">$outputpath/results/novel_junction_details.tab");
    print NOVELDETAILS "subgenome\tleader_end\tTRS_start\tACGAAC\tATG_postion\tknown_ATG\t20_leader_seq\tTRS_seq\tfirst_orf_aa\n";
    open (KNOWNJUNCATIONSDETAILS, ">$outputpath/results/known_junction_details.tab");
    print KNOWNJUNCATIONSDETAILS "subgenome\tpeak_leader_end\tpeak_TRS_start\tACGAAC\tATG_postion\t20_leader_seq\tTRS_seq\tfirst_orf_aa\n";

    open (NOVEL, ">$outputpath/results/novel_junction.tab");
    print NOVEL "subgenome\tleader_end\tTRS_start\tnb_count\tnormalized_count\n";
    
    open (KNOWNJUNCATIONS, ">$outputpath/results/known_junction.tab");
    print KNOWNJUNCATIONS "subgenome\tref_leader_end\tpeak_leader_end\tref_TRS_start\tpeak_TRS_start\tpeak_count\tpeak_normalized_count\tcluster_count\tcluster_normalized_count\n";
    
    open (PRIMERMOSTUSED, ">$outputpath/results/primers_used.tab");
    
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
            if ($eachtabs[4]<$cutoffpos && $eachtabs[20]>$cutoff) {
                $hashtab{$eachtabs[0]}="$eachtabs[0]\t$eachtabs[4]\t$eachtabs[5]\t$eachtabs[20]";
                push (@allindexcollection,"$eachtabs[0]");
            }
            #if ($eachtabs[4]<266) {
                $hashtabr{$eachtabs[0]}="$eachtabs[0]\t$eachtabs[4]\t$eachtabs[5]\t$eachtabs[20]";
            #}
        }
    }
    close TAB;
    
    open(BAMEXTRACTR, "$outputpath/alignment_output/alignment.sorted.splice.tab");
    my %hashassp; my %hashasspleftright; my %hashassfa; my %hashillumina; my %hashasspleftright1; my %hashasspleftright2; #%hashassp for right site; %hashasspleftright for left and right site; %hashassfa for fasta extaction;
    while (<BAMEXTRACTR>) {
        my @eachassps=split(/\t/);
        if ($eachassps[1] eq 1) {
            if ($options{'mode'} eq "illumia" && !exists $options{'bam'}) {
                push (@{$hashassp{"$eachassps[3]\t$eachassps[4]"}}, "$eachassps[5]");
                my @eachleftright1=($eachassps[0],$eachassps[2],$eachassps[5],$eachassps[9],$eachassps[11],$eachassps[14]);
                push (@{$hashasspleftright1{"$eachassps[3]\t$eachassps[4]"}}, [@eachleftright1]);
                my @eachleftright2=($eachassps[0],$eachassps[2],$eachassps[5],$eachassps[9],$eachassps[11],$eachassps[14]);
                push (@{$hashasspleftright2{"$eachassps[12]\t$eachassps[13]"}}, [@eachleftright2]);
                
                push (@{$hashassfa{"$eachassps[3]\t$eachassps[4]"}}, "$eachassps[0]\t$eachassps[2]\t$eachassps[5]\t$eachassps[8]");
                
                my @spliteachassps0=split(/\_\_/,$eachassps[0]);
                push (@{$hashillumina{"$eachassps[3]\t$eachassps[4]"}}, "$spliteachassps0[1]");
            }else {
                push (@{$hashassp{"$eachassps[3]\t$eachassps[4]"}}, "$eachassps[5]");
                my @eachleftright=($eachassps[0],$eachassps[2],$eachassps[5]);
                push (@{$hashasspleftright{"$eachassps[3]\t$eachassps[4]"}}, [@eachleftright]);
                push (@{$hashassfa{"$eachassps[3]\t$eachassps[4]"}}, "$eachassps[0]\t$eachassps[2]\t$eachassps[5]\t$eachassps[8]");
            }
        }
    }
    close BAMEXTRACTR;
    
    my @indexcollections;

    foreach my $junction (@junctions) {
        chomp $junction;
        my @junctionsites=split(/\t/,$junction);
        $junctionhash{$junctionsites[3]}=$junctionsites[0];
        my %hashcluster=();
        my %hashleftrightprimerscluster=();
        my %hashleftrightprimerspeak=();
        
        my @clustercounts=(); #for the normal counts collection
        
        my @clusterployAcounts1=(); #for the ployA >= 1 counts collection
        my @clusterployAcounts5=(); #for the ployA >= 5 counts collection
        
        my @clusterleftprimerids=(); #for the left primer ids collection
        my @clusterrightprimerids=(); #for the right primer ids collection
        
        my %hashclusterfirstleftrightprimerids=();
        
        #my @clusterATGsite=(); #for the ATGids collection
        
        #####cluster count of junction
        open(TAB, "$outputpath/alignment_output/alignment_portcullis_out.junctions.tab");
        while (<TAB>) {
            my @eachtabs=split(/\t/);
            unless (/^index\trefid/ or /^\n/){
                my $adjectlederup=$junctionsites[1]+$leaderadjectnumber;
                my $adjectlederdown=$junctionsites[1]-$leaderadjectnumber;
                my $adjectorfup=$junctionsites[3];
                my $adjectorfdown=$junctionsites[2]-$orfadjectnumber;
                if ($eachtabs[4] < $adjectlederup && $eachtabs[4] > $adjectlederdown && $eachtabs[5] < $adjectorfup && $eachtabs[5] > $adjectorfdown) {
                    push (@indexcollections,$eachtabs[0]);
                    push (@clustercounts,$eachtabs[20]); #cluster uniqe junction count.
                    $hashcluster{$eachtabs[0]}=$eachtabs[20];
                    
                    if ($options{'Rtch'} eq "RNA" && !exists $options{'bam'}) {
                        my $countployAcounts1=grep($_>=$polyapostion, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        push (@clusterployAcounts1,$countployAcounts1); # to get ployA >= 1 by the sart and end postions in the cluster.
                        my $countployAcounts5=grep($_>=$polyapostion5, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        push (@clusterployAcounts5,$countployAcounts5); # to get ployA >=5 by the sart and end postions in the cluster.
                        #my $clusterATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        #push (@clusterATGsite,$clusterATGsiteid);
                    }
                    
                    if ($options{'Rtch'} eq "cDNA" && !exists $options{'bam'}) {
                        #my $countployAcounts1=grep($_>=$polyapostion, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        #push (@clusterployAcounts1,$countployAcounts1); # to get ployA >= 1 by the sart and end postions in the cluster.
                        #my $countployAcounts5=grep($_>=$polyapostion5, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        #push (@clusterployAcounts5,$countployAcounts5); # to get ployA >=5 by the sart and end postions in the cluster.
                        #my $clusterATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabs[4]\t$eachtabs[5]"}});
                        #push (@clusterATGsite,$clusterATGsiteid);
                        
                        open(PRIMERBED, "$options{'primer_bed'}");
                        while (<PRIMERBED>) {
                            chomp;
                            my @eachprimerbed=split(/\t/);
                            if ($eachprimerbed[1] < $eachtabs[4] && $eachprimerbed[3]=~/_LEFT/) {
                                my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachtabs[4]\t$eachtabs[5]"}});
                                for (my $n=0; $n<$#leftprimerids+1; $n++) {
                                    push (@clusterleftprimerids,$leftprimerids[$n][0]);
                                    push (@{$hashclusterfirstleftrightprimerids{$leftprimerids[$n][0]}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                                }
                            }
                                                        
                            if ($eachprimerbed[2] > $eachtabs[5] && $eachprimerbed[3]=~/_RIGHT/) {
                                my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachtabs[4]\t$eachtabs[5]"}});
                                for (my $n=0; $n<$#rightprimerids+1; $n++) {
                                    push (@clusterrightprimerids,$rightprimerids[$n][0]);
                                    push (@{$hashclusterfirstleftrightprimerids{$rightprimerids[$n][0]}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
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
                            if ($eachprimerbed[1] < $eachtabs[4] && $eachprimerbed[3]=~/_LEFT/) {   
                                my @leftprimerids1=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright1{"$eachtabs[4]\t$eachtabs[5]"}});
                                my @leftprimerids2=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright2{"$eachtabs[4]\t$eachtabs[5]"}});
                                for (my $n=0; $n<$#leftprimerids1+1; $n++) {
                                    push (@clusterleftprimerids,"$leftprimerids1[$n][0]\_$leftprimerids1[$n][3]");
                                    push (@{$hashclusterfirstleftrightprimerids{"$leftprimerids1[$n][0]\_$leftprimerids1[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                                }
                                for (my $n=0; $n<$#leftprimerids2+1; $n++) {
                                    push (@clusterleftprimerids,"$leftprimerids2[$n][0]\_$leftprimerids2[$n][3]");
                                    push (@{$hashclusterfirstleftrightprimerids{"$leftprimerids2[$n][0]\_$leftprimerids2[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                                }
                            }
                                                        
                            if ($eachprimerbed[2] > $eachtabs[5] && $eachprimerbed[3]=~/_RIGHT/) {
                                my @rightprimerids1=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[5]), @{$hashasspleftright1{"$eachtabs[4]\t$eachtabs[5]"}});
                                my @rightprimerids2=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[5]), @{$hashasspleftright2{"$eachtabs[4]\t$eachtabs[5]"}});
                                for (my $n=0; $n<$#rightprimerids1+1; $n++) {
                                    push (@clusterrightprimerids,"$rightprimerids1[$n][0]\_$rightprimerids1[$n][3]");
                                    push (@{$hashclusterfirstleftrightprimerids{"$rightprimerids1[$n][0]\_$rightprimerids1[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
                                }
                                for (my $n=0; $n<$#rightprimerids2+1; $n++) {
                                    push (@clusterrightprimerids,"$rightprimerids2[$n][0]\_$rightprimerids2[$n][3]");
                                    push (@{$hashclusterfirstleftrightprimerids{"$rightprimerids2[$n][0]\_$rightprimerids2[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
                                }
                            }
                        }
                        close PRIMERBED;
                    }
                }
            }
        }
        close TAB;
        
        my @clusters= (sort {$hashcluster{$a} <=> $hashcluster{$b}} keys %hashcluster);
        
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
        my $clusterleftrightprimeridsunioncount=$#clusterleftrightprimeridsunion+1;
        
        my @clusterfirstleftrightprimeridsclustesrpro=();
        foreach (@clusterleftrightprimeridsintersection) {
            push (@clusterfirstleftrightprimeridsclustesrpro,join("\_",uniq @{$hashclusterfirstleftrightprimerids{$_}}));
        }
        my %countsfirstleftrightprimerids = ();
        for (@clusterfirstleftrightprimeridsclustesrpro) {
            $countsfirstleftrightprimerids{$_}++;
        }
        my @clusterfirstleftrightprimeridsclusters= (sort {$countsfirstleftrightprimerids{$a} <=> $countsfirstleftrightprimerids{$b}} keys %countsfirstleftrightprimerids);
        
        my $clusterfirstleftrightprimeridsintersectioncount=0;
        if ($#clusterfirstleftrightprimeridsclusters == -1) {
            $clusterfirstleftrightprimeridsintersectioncount=0;
        }else{
            $clusterfirstleftrightprimeridsintersectioncount=$countsfirstleftrightprimerids{$clusterfirstleftrightprimeridsclusters[-1]};
            print PRIMERMOSTUSED "$junctionsites[0]\t$clusterfirstleftrightprimeridsclusters[-1]\t$clusterfirstleftrightprimeridsintersectioncount\n";
        }

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
                my $countployAcounts1=grep($_>=$polyapostion, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                my $countployAcounts5=grep($_>=$polyapostion5, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
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
                #my $covsumclusterployAcounts1=sprintf("%.2f", $sumclusterployAcounts1/$mappedcoverage[0]*1000000);
                #my $covsumclusterployAcounts5=sprintf("%.2f", $sumclusterployAcounts5/$mappedcoverage[0]*1000000);
                
                #my $countATGsiteid=grep($_>$adjectorfup+1, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                #my $countployAcounts1=grep($_>=$polyapostion, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                #my $countployAcounts5=grep($_>=$polyapostion5, @{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                #my $covcountATGsiteid=sprintf("%.2f", $countATGsiteid/$mappedcoverage[0]*1000000);
                #my $covcountployAcounts1=sprintf("%.2f", $countployAcounts1/$mappedcoverage[0]*1000000);
                #my $covcountployAcounts5=sprintf("%.2f", $countployAcounts5/$mappedcoverage[0]*1000000);
                
                my @clusterleftprimeridspeak=();
                my @clusterrightprimeridspeak=();
                my %hashclusterfirstleftrightprimeridspeak=();
                
                open(PRIMERBED, "$options{'primer_bed'}");
                while (<PRIMERBED>) {
                    chomp;
                    my @eachprimerbed=split(/\t/);
                    if ($eachprimerbed[1] < $eachtabsagain[1] && $eachprimerbed[3]=~/_LEFT/) {
                        my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        for (my $n=0; $n<$#leftprimerids+1; $n++) {
                            push (@clusterleftprimeridspeak,$leftprimerids[$n][0]);
                            push (@{$hashclusterfirstleftrightprimeridspeak{$leftprimerids[$n][0]}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                        }
                    }
           
                    if ($eachprimerbed[2] > $eachtabsagain[2] && $eachprimerbed[3]=~/_RIGHT/) {   
                        my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        for (my $n=0; $n<$#rightprimerids+1; $n++) {
                            push (@clusterrightprimeridspeak,$rightprimerids[$n][0]);
                            push (@{$hashclusterfirstleftrightprimeridspeak{$rightprimerids[$n][0]}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
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
                
                
                my @clusterfirstleftrightprimeridsclustesrpropeak=();
                foreach (@clusterleftrightprimeridsintersectionpeak) {
                    push (@clusterfirstleftrightprimeridsclustesrpropeak,join("\_",uniq @{$hashclusterfirstleftrightprimeridspeak{$_}}));
                }
                my %countsfirstleftrightprimeridspeak = ();
                for (@clusterfirstleftrightprimeridsclustesrpropeak) {
                    $countsfirstleftrightprimeridspeak{$_}++;
                }
                my @clusterfirstleftrightprimeridsclusterspeak= (sort {$countsfirstleftrightprimeridspeak{$a} <=> $countsfirstleftrightprimeridspeak{$b}} keys %countsfirstleftrightprimeridspeak);
                my $clusterfirstleftrightprimeridsintersectioncountpeak=0;
                if ($#clusterfirstleftrightprimeridsclusterspeak == -1) {
                    $clusterfirstleftrightprimeridsintersectioncountpeak=0;
                }else{
                    $clusterfirstleftrightprimeridsintersectioncountpeak=$countsfirstleftrightprimeridspeak{$clusterfirstleftrightprimeridsclusterspeak[-1]};
                    print PRIMERMOSTUSED "$junctionsites[0]\t$clusterfirstleftrightprimeridsclusterspeak[-1]\t$clusterfirstleftrightprimeridsintersectioncountpeak\n";
                }
                
                my $covsumclusterleftprimercountspeak=sprintf("%.2f", $sumclusterleftprimercountspeak/$mappedcoverage[0]*1000000);
                my $covsumclusterrightprimercountspeak=sprintf("%.2f", $sumclusterrightprimercountspeak/$mappedcoverage[0]*1000000);
                my $covclusterleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]*1000000);
                my $covclusterfirstleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterfirstleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]*1000000);
                
                my $covsumclusterleftprimercounts=sprintf("%.2f", $sumclusterleftprimercounts/$mappedcoverage[0]*1000000);
                my $covsumclusterrightprimercounts=sprintf("%.2f", $sumclusterrightprimercounts/$mappedcoverage[0]*1000000);
                my $covclusterleftrightprimeridsintersectioncount=sprintf("%.2f", $clusterleftrightprimeridsintersectioncount/$mappedcoverage[0]*1000000);
                my $covclusterfirstleftrightprimeridsintersectioncount=sprintf("%.2f", $clusterfirstleftrightprimeridsintersectioncount/$mappedcoverage[0]*1000000);
                
                #my $testlength=$#{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}}+1;
                &extractdetails ($junctionsites[0],$eachtabsagain[1],$eachtabsagain[2]);
                #print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t$eachtabsagain[1]\t$junctionsites[2]\t$eachtabsagain[2]\t$eachtabsagain[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$countployAcounts1,$countployAcounts5\)\t$covratiopeack\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covcountployAcounts1,$covcountployAcounts5\)\t$sumclustercounts\($sumclusterleftprimercounts,$sumclusterrightprimercounts,$clusterleftrightprimeridsintersectioncount,$sumclusterployAcounts1,$sumclusterployAcounts5\)\t$covratiocluster\($covsumclusterleftprimercounts,$covsumclusterrightprimercounts,$covclusterleftrightprimeridsintersectioncount,$covsumclusterployAcounts1,$covsumclusterployAcounts5\)\n";
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t$eachtabsagain[1]\t$junctionsites[2]\t$eachtabsagain[2]\t$eachtabsagain[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$clusterfirstleftrightprimeridsintersectioncountpeak\)\t$covratiopeack\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covclusterfirstleftrightprimeridsintersectioncountpeak\)\t$sumclustercounts\($sumclusterleftprimercounts,$sumclusterrightprimercounts,$clusterleftrightprimeridsintersectioncount,$clusterfirstleftrightprimeridsintersectioncount\)\t$covratiocluster\($covsumclusterleftprimercounts,$covsumclusterrightprimercounts,$covclusterleftrightprimeridsintersectioncount,$covclusterfirstleftrightprimeridsintersectioncount\)\n";
            }
        }elsif ($options{'mode'} eq "illumia" && !exists $options{'bam'}) {
            if ($#clusters == -1) {
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t0\t$junctionsites[2]\t0\t0\t0\t0\t0\n";
            }else{
                my @eachtabsagain=split(/\t/,$hashtabr{$clusters[-1]});

                my @clusterleftprimeridspeak=();
                my @clusterrightprimeridspeak=();
                my %hashclusterfirstleftrightprimeridspeak=();
                
                open(PRIMERBED, "$options{'primer_bed'}");
                while (<PRIMERBED>) {
                    chomp;
                    my @eachprimerbed=split(/\t/);
                    if ($eachprimerbed[1] < $eachtabsagain[1] && $eachprimerbed[3]=~/_LEFT/) {   
                        my @leftprimerids1=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright1{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        my @leftprimerids2=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright2{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        for (my $n=0; $n<$#leftprimerids1+1; $n++) {
                            push (@clusterleftprimeridspeak,"$leftprimerids1[$n][0]\_$leftprimerids1[$n][3]");
                            push (@{$hashclusterfirstleftrightprimeridspeak{"$leftprimerids1[$n][0]\_$leftprimerids1[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                        }
                        for (my $n=0; $n<$#leftprimerids2+1; $n++) {
                            push (@clusterleftprimeridspeak,"$leftprimerids2[$n][0]\_$leftprimerids2[$n][3]");
                            push (@{$hashclusterfirstleftrightprimeridspeak{"$leftprimerids2[$n][0]\_$leftprimerids2[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                        }
                    }
                                                        
                    if ($eachprimerbed[2] > $eachtabsagain[2] && $eachprimerbed[3]=~/_RIGHT/) {
                        my @rightprimerids1=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[5]), @{$hashasspleftright1{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        my @rightprimerids2=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[5]), @{$hashasspleftright2{"$eachtabsagain[1]\t$eachtabsagain[2]"}});
                        for (my $n=0; $n<$#rightprimerids1+1; $n++) {
                            push (@clusterrightprimeridspeak,"$rightprimerids1[$n][0]\_$rightprimerids1[$n][3]");
                            push (@{$hashclusterfirstleftrightprimeridspeak{"$rightprimerids1[$n][0]\_$rightprimerids1[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
                        }
                        for (my $n=0; $n<$#rightprimerids2+1; $n++) {
                            push (@clusterrightprimeridspeak,"$rightprimerids2[$n][0]\_$rightprimerids2[$n][3]");
                            push (@{$hashclusterfirstleftrightprimeridspeak{"$rightprimerids2[$n][0]\_$rightprimerids2[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
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
                my $clusterleftrightprimeridsunioncountpeak=$#clusterleftrightprimeridsunionpeak+1;
                
                my @clusterfirstleftrightprimeridsclustesrpropeak=();
                foreach (@clusterleftrightprimeridsintersectionpeak) {
                    push (@clusterfirstleftrightprimeridsclustesrpropeak,join("\_",uniq @{$hashclusterfirstleftrightprimeridspeak{$_}}));
                }
                my %countsfirstleftrightprimeridspeak = ();
                for (@clusterfirstleftrightprimeridsclustesrpropeak) {
                    $countsfirstleftrightprimeridspeak{$_}++;
                }
                my @clusterfirstleftrightprimeridsclusterspeak= (sort {$countsfirstleftrightprimeridspeak{$a} <=> $countsfirstleftrightprimeridspeak{$b}} keys %countsfirstleftrightprimeridspeak);
                my $clusterfirstleftrightprimeridsintersectioncountpeak=0;
                if ($#clusterfirstleftrightprimeridsclusterspeak == -1) {
                    $clusterfirstleftrightprimeridsintersectioncountpeak=0;
                }else{
                    $clusterfirstleftrightprimeridsintersectioncountpeak=$countsfirstleftrightprimeridspeak{$clusterfirstleftrightprimeridsclusterspeak[-1]};
                    print PRIMERMOSTUSED "$junctionsites[0]\t$clusterfirstleftrightprimeridsclusterspeak[-1]\t$clusterfirstleftrightprimeridsintersectioncountpeak\n";
                }
                
                my $covsumclusterleftprimercountspeak=sprintf("%.2f", $sumclusterleftprimercountspeak/$mappedcoverage[0]/2*1000000);
                my $covsumclusterrightprimercountspeak=sprintf("%.2f", $sumclusterrightprimercountspeak/$mappedcoverage[0]/2*1000000);
                my $covclusterleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]/2*1000000);
                my $covclusterleftrightprimeridsunioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsunioncountpeak/$mappedcoverage[0]/2*1000000);
                my $covclusterfirstleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterfirstleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]/2*1000000);
                
                my $covsumclusterleftprimercounts=sprintf("%.2f", $sumclusterleftprimercounts/$mappedcoverage[0]/2*1000000);
                my $covsumclusterrightprimercounts=sprintf("%.2f", $sumclusterrightprimercounts/$mappedcoverage[0]/2*1000000);
                my $covclusterleftrightprimeridsintersectioncount=sprintf("%.2f", $clusterleftrightprimeridsintersectioncount/$mappedcoverage[0]/2*1000000);
                my $covclusterleftrightprimeridsunioncount=sprintf("%.2f", $clusterleftrightprimeridsunioncount/$mappedcoverage[0]/2*1000000);
                my $covclusterfirstleftrightprimeridsintersectioncount=sprintf("%.2f", $clusterfirstleftrightprimeridsintersectioncount/$mappedcoverage[0]/2*1000000);
            

                ########my $numberfirstpairedjucntioncluster
                
                #my $testlength=$#{$hashassp{"$eachtabsagain[1]\t$eachtabsagain[2]"}}+1;
                &extractdetails ($junctionsites[0],$eachtabsagain[1],$eachtabsagain[2]);
                #print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t$eachtabsagain[1]\t$junctionsites[2]\t$eachtabsagain[2]\t$eachtabsagain[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$numberpairedjucntionpeak\)\t$covratiopeack\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covnumberpairedjucntionpeak\)\t$sumclustercounts\($sumclusterleftprimercounts,$sumclusterrightprimercounts,$clusterleftrightprimeridsintersectioncount,$numberpairedjucntioncluster\)\t$covratiocluster\($covsumclusterleftprimercounts,$covsumclusterrightprimercounts,$covclusterleftrightprimeridsintersectioncount,$covnumberpairedjucntioncluster\)\n";
                print KNOWNJUNCATIONS "$junctionsites[0]\t$junctionsites[1]\t$eachtabsagain[1]\t$junctionsites[2]\t$eachtabsagain[2]\t$clusterleftrightprimeridsunioncountpeak\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$clusterfirstleftrightprimeridsintersectioncountpeak\)\t$covclusterleftrightprimeridsunioncountpeak\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covclusterfirstleftrightprimeridsintersectioncountpeak\)\t$clusterleftrightprimeridsunioncount\($sumclusterleftprimercounts,$sumclusterrightprimercounts,$clusterleftrightprimeridsintersectioncount,$clusterfirstleftrightprimeridsintersectioncount\)\t$covclusterleftrightprimeridsunioncount\($covsumclusterleftprimercounts,$covsumclusterrightprimercounts,$covclusterleftrightprimeridsintersectioncount,$covclusterfirstleftrightprimeridsintersectioncount\)\n";
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

    ########novel subgenomes
    my $indexnovel = List::Compare->new(\@indexcollections, \@allindexcollection);
    
    my @indexnovelcollections = $indexnovel->get_complement;
    
    if ($options{'Rtch'} eq "RNA" && !exists $options{'bam'}) {
        my $n=0;
        foreach my $indexnovelcollection(@indexnovelcollections) {
            $n++;
            my @eachnvelintable=split(/\t/,"$hashtab{$indexnovelcollection}");
            my $countployAcounts1novel=grep($_>=$polyapostion, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            my $countployAcounts5novel=grep($_>=$polyapostion5, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
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
            #my $countployAcounts1novel=grep($_>=$polyapostion, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            #my $countployAcounts5novel=grep($_>=$polyapostion5, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            #my $countATGsiteid=grep($_>$atgporstionnovel+1, @{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
            my $covnoveluniqe=sprintf("%.2f", $eachnvelintable[3]/$mappedcoverage[0]*1000000);
            #my $covcountployAcounts1novel=sprintf("%.2f", $countployAcounts1novel/$mappedcoverage[0]*1000000);
            #my $covcountployAcounts5novel=sprintf("%.2f", $countployAcounts5novel/$mappedcoverage[0]*1000000);
            #my $covcountATGsiteid=sprintf("%.2f", $countATGsiteid/$mappedcoverage[0]*1000000);
            
            my @clusterleftprimeridspeak=();
            my @clusterrightprimeridspeak=();
            my %hashclusterfirstleftrightprimeridspeak=();
            open(PRIMERBED, "$options{'primer_bed'}");
            while (<PRIMERBED>) {
                chomp;
                my @eachprimerbed=split(/\t/);
                if ($eachprimerbed[1] < $eachnvelintable[1] && $eachprimerbed[3]=~/_LEFT/) {
                    my @leftprimerids=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    for (my $n=0; $n<$#leftprimerids+1; $n++) {
                        push (@clusterleftprimeridspeak,$leftprimerids[$n][0]);
                        push (@{$hashclusterfirstleftrightprimeridspeak{$leftprimerids[$n][0]}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                    }
                }
                                
                if ($eachprimerbed[2] > $eachnvelintable[2] && $eachprimerbed[3]=~/_RIGHT/) {     
                    my @rightprimerids=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[2]), @{$hashasspleftright{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    for (my $n=0; $n<$#rightprimerids+1; $n++) {
                        push (@clusterrightprimeridspeak,$rightprimerids[$n][0]);
                        push (@{$hashclusterfirstleftrightprimeridspeak{$rightprimerids[$n][0]}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
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
                        
            my @clusterfirstleftrightprimeridsclustesrpropeak=();
            foreach (@clusterleftrightprimeridsintersectionpeak) {
                push (@clusterfirstleftrightprimeridsclustesrpropeak,join("\_",uniq @{$hashclusterfirstleftrightprimeridspeak{$_}}));
            }
            my %countsfirstleftrightprimeridspeak = ();
            for (@clusterfirstleftrightprimeridsclustesrpropeak) {
                $countsfirstleftrightprimeridspeak{$_}++;
            }
            my @clusterfirstleftrightprimeridsclusterspeak= (sort {$countsfirstleftrightprimeridspeak{$a} <=> $countsfirstleftrightprimeridspeak{$b}} keys %countsfirstleftrightprimeridspeak);
            my $clusterfirstleftrightprimeridsintersectioncountpeak=0;
            if ($#clusterfirstleftrightprimeridsclusterspeak == -1) {
                $clusterfirstleftrightprimeridsintersectioncountpeak=0;
            }else{
                $clusterfirstleftrightprimeridsintersectioncountpeak=$countsfirstleftrightprimeridspeak{$clusterfirstleftrightprimeridsclusterspeak[-1]};
                print PRIMERMOSTUSED "$n\_novel\t$clusterfirstleftrightprimeridsclusterspeak[-1]\t$clusterfirstleftrightprimeridsintersectioncountpeak\n";
            }
                
            my $covsumclusterleftprimercountspeak=sprintf("%.2f", $sumclusterleftprimercountspeak/$mappedcoverage[0]*1000000);
            my $covsumclusterrightprimercountspeak=sprintf("%.2f", $sumclusterrightprimercountspeak/$mappedcoverage[0]*1000000);
            my $covclusterleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]*1000000);
            my $covclusterfirstleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterfirstleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]*1000000);
            
            #my $testlength=$#{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}}+1;
            
            &extractdetailsnovel($n,$eachnvelintable[1],$eachnvelintable[2]);
            #print NOVEL "$n\t$eachnvelintable[1]\t$eachnvelintable[2]\t$eachnvelintable[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$countployAcounts1novel,$countployAcounts5novel\)\t$covnoveluniqe\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covcountployAcounts1novel,$covcountployAcounts5novel\)\n";
            print NOVEL "$n\t$eachnvelintable[1]\t$eachnvelintable[2]\t$eachnvelintable[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$clusterfirstleftrightprimeridsintersectioncountpeak\)\t$covnoveluniqe\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covclusterfirstleftrightprimeridsintersectioncountpeak\)\n";
        }
    }elsif ($options{'mode'} eq "illumia" && !exists $options{'bam'}) {
        my $n=0;
        foreach my $indexnovelcollection(@indexnovelcollections) {
            $n++;
            my @eachnvelintable=split(/\t/,"$hashtab{$indexnovelcollection}");
          
            my @clusterleftprimeridspeak=();
            my @clusterrightprimeridspeak=();
            my %hashclusterfirstleftrightprimeridspeak=();
                
            open(PRIMERBED, "$options{'primer_bed'}");
            while (<PRIMERBED>) {
                chomp;
                my @eachprimerbed=split(/\t/);
                if ($eachprimerbed[1] < $eachnvelintable[1] && $eachprimerbed[3]=~/_LEFT/) {   
                    my @leftprimerids1=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright1{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    my @leftprimerids2=grep(grep($_>=$eachprimerbed[1] && $_<=$eachprimerbed[2]-10, @$_[1]), @{$hashasspleftright2{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    for (my $n=0; $n<$#leftprimerids1+1; $n++) {
                        push (@clusterleftprimeridspeak,"$leftprimerids1[$n][0]\_$leftprimerids1[$n][3]");
                        push (@{$hashclusterfirstleftrightprimeridspeak{"$leftprimerids1[$n][0]\_$leftprimerids1[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                    }
                    for (my $n=0; $n<$#leftprimerids2+1; $n++) {
                        push (@clusterleftprimeridspeak,"$leftprimerids2[$n][0]\_$leftprimerids2[$n][3]");
                        push (@{$hashclusterfirstleftrightprimeridspeak{"$leftprimerids2[$n][0]\_$leftprimerids2[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_left");
                    }
                }
                                                        
                if ($eachprimerbed[2] > $eachnvelintable[2] && $eachprimerbed[3]=~/_RIGHT/) {
                    my @rightprimerids1=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[5]), @{$hashasspleftright1{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    my @rightprimerids2=grep(grep($_>=$eachprimerbed[1]+10 && $_<=$eachprimerbed[2], @$_[5]), @{$hashasspleftright2{"$eachnvelintable[1]\t$eachnvelintable[2]"}});
                    for (my $n=0; $n<$#rightprimerids1+1; $n++) {
                        push (@clusterrightprimeridspeak,"$rightprimerids1[$n][0]\_$rightprimerids1[$n][3]");
                        push (@{$hashclusterfirstleftrightprimeridspeak{"$rightprimerids1[$n][0]\_$rightprimerids1[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
                    }
                    for (my $n=0; $n<$#rightprimerids2+1; $n++) {
                        push (@clusterrightprimeridspeak,"$rightprimerids2[$n][0]\_$rightprimerids2[$n][3]");
                        push (@{$hashclusterfirstleftrightprimeridspeak{"$rightprimerids2[$n][0]\_$rightprimerids2[$n][3]"}},"$eachprimerbed[1]_$eachprimerbed[2]_right");
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
            my $clusterleftrightprimeridsunioncountpeak=$#clusterleftrightprimeridsunionpeak+1;
            
            my @clusterfirstleftrightprimeridsclustesrpropeak=();
            foreach (@clusterleftrightprimeridsintersectionpeak) {
                push (@clusterfirstleftrightprimeridsclustesrpropeak,join("\_",uniq @{$hashclusterfirstleftrightprimeridspeak{$_}}));
            }
            my %countsfirstleftrightprimeridspeak = ();
            for (@clusterfirstleftrightprimeridsclustesrpropeak) {
                $countsfirstleftrightprimeridspeak{$_}++;
            }
            my @clusterfirstleftrightprimeridsclusterspeak= (sort {$countsfirstleftrightprimeridspeak{$a} <=> $countsfirstleftrightprimeridspeak{$b}} keys %countsfirstleftrightprimeridspeak);
            my $clusterfirstleftrightprimeridsintersectioncountpeak=0;
            if ($#clusterfirstleftrightprimeridsclusterspeak == -1) {
                $clusterfirstleftrightprimeridsintersectioncountpeak=0;
            }else{
                $clusterfirstleftrightprimeridsintersectioncountpeak=$countsfirstleftrightprimeridspeak{$clusterfirstleftrightprimeridsclusterspeak[-1]};
                print PRIMERMOSTUSED "$n\_novel\t$clusterfirstleftrightprimeridsclusterspeak[-1]\t$clusterfirstleftrightprimeridsintersectioncountpeak\n";
            }
                
            my $covsumclusterleftprimercountspeak=sprintf("%.2f", $sumclusterleftprimercountspeak/$mappedcoverage[0]/2*1000000);
            my $covsumclusterrightprimercountspeak=sprintf("%.2f", $sumclusterrightprimercountspeak/$mappedcoverage[0]/2*1000000);
            my $covclusterleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]/2*1000000);
            my $covclusterleftrightprimeridsunioncountpeak=sprintf("%.2f", $clusterleftrightprimeridsunioncountpeak/$mappedcoverage[0]/2*1000000);
            my $covclusterfirstleftrightprimeridsintersectioncountpeak=sprintf("%.2f", $clusterfirstleftrightprimeridsintersectioncountpeak/$mappedcoverage[0]/2*1000000);
   
            #my $testlength=$#{$hashassp{"$eachnvelintable[1]\t$eachnvelintable[2]"}}+1;
            &extractdetailsnovel($n,$eachnvelintable[1],$eachnvelintable[2]);
            #print NOVEL "$n\t$eachnvelintable[1]\t$eachnvelintable[2]\t$eachnvelintable[3]\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$numberpairedjucntionpeak\)\t$covnoveluniqe\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covnumberpairedjucntionpeak\)\n";
            print NOVEL "$n\t$eachnvelintable[1]\t$eachnvelintable[2]\t$clusterleftrightprimeridsunioncountpeak\($sumclusterleftprimercountspeak,$sumclusterrightprimercountspeak,$clusterleftrightprimeridsintersectioncountpeak,$clusterfirstleftrightprimeridsintersectioncountpeak\)\t$covclusterleftrightprimeridsunioncountpeak\($covsumclusterleftprimercountspeak,$covsumclusterrightprimercountspeak,$covclusterleftrightprimeridsintersectioncountpeak,$covclusterfirstleftrightprimeridsintersectioncountpeak\)\n";
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
        print KNOWNJUNCATIONS "The numbers in the bracket are (reads with left primers, reads with right primers, reads with both primers, reads with most used primers)\.\n";
        print KNOWNJUNCATIONS "Normalized count=(Read count/Total number of read mapped on reference genome)*1000000\.\n";
        print KNOWNJUNCATIONS "Total number of read mapped on reference genome is $mappedcoverage[0], excluding the mapped reads not primary alignment and supplementary alignment\.\n";
        print NOVEL "The numbers in the bracket are (reads with left primers, reads with right primers, reads with both primers, reads with most used primers)\.\n";
        print NOVEL "Normalized count=(Read count/Total number of read mapped on reference genome)*1000000\.\n";
        print NOVEL "Total number of read mapped on reference genome doesn't include the mapped reads not primary alignment and supplementary alignment\.\n";
    }elsif ($options{'mode'} eq "illumia" && !exists $options{'bam'}) {
        my $mappedcoveragenew=sprintf("%.0f",$mappedcoverage[0]/2);
        print KNOWNJUNCATIONS "The numbers in the bracket are (read pairs with left primers, read pairs with right primers, read pairs with both primers, read pairs with most used primers)\.\n";
        print KNOWNJUNCATIONS "Normalized count=(Read pair count/Total number of read pair mapped on reference genome)*1000000\.\n";
        print KNOWNJUNCATIONS "Total number of read pair mapped on reference genome is $mappedcoveragenew, excluding the mapped reads unpaired, not primary alignment and supplementary alignment\.\n";
        print NOVEL "The numbers in the bracket are (read pairs with left primers, read pairs with right primers, read pairs with both primers, read pairs with most used primers)\.\n";
        print NOVEL "Normalized count=(Read pair count/Total number of read pair mapped on reference genome)*1000000\.\n";
        print NOVEL "Total number of read pair mapped on reference is $mappedcoveragenew, excluding the mapped reads unpaired, not primary alignment and supplementary alignment\.\n";
    }else{
        print KNOWNJUNCATIONS "Total number of read mapped on reference genome is $mappedcoverage[0], excluding the mapped reads unpaired, not primary alignment and supplementary alignment\.\n";
        print NOVEL "Total number of read mapped on reference is $mappedcoverage[0], excluding the mapped reads unpaired, not primary alignment and supplementary alignment\.\n";
    }
    close NOVEL; close KNOWNJUNCATIONS;
    close KNOWNJUNCATIONSDETAILS; close NOVELDETAILS; close PRIMERMOSTUSED;
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

            my $prot_obj =$seq->trunc($each4-1,$genomelength)->translate(-orf => 1, -start => "atg");
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
        if (defined pos($TRSseqall)) {
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

            my $prot_obj =$seq->trunc($each2-1,$genomelength)->translate(-orf => 1, -start => "atg");
            print NOVELDETAILS "$atgporstion\t";
            print NOVELDETAILS "$known_ATG\t";
            print NOVELDETAILS "$leaderseq\t";
            print NOVELDETAILS "$TRSseq\t";
            print NOVELDETAILS $prot_obj->seq, "\n";
        }else {
            print NOVELDETAILS "\-\t";
            print NOVELDETAILS "\-\t";
            print NOVELDETAILS "\-\t";
            print NOVELDETAILS "\-\n";
            print NOVELDETAILS "\-\n";
        }
    }
}
