#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my %options;
my @standard_options =("help|h!",
                       "genome=s",
                       "bed=s",
                       "o=s",
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
  perl making_ref.pl -genome genome.fasta -bed junction.bed -o reference_folder
   
  -genome           input fasta file of a coronavirus reference genome.
  -bed              a four columns bed file including subgenomic name, position of the end of leader,
                    position of the start of TRS and position of the end of TRS.
  -o                output path, \"./\" by default.
                    
  -h/-help          Produce help message.\n\n";
    exit;
}

if (!exists $options{'genome'}) {
    print "please provide a coronavirus reference genome in fasta format";
}
if (!exists $options{'bed'}) {
    print "please provide a four columns bed file including subgenomic name, position of the end of leader, position of the start of TRS and position of the end of TRS.";
}
    
my $outputpath;
if (exists $options{'o'}) {
   $outputpath=$options{'o'};  
}else{
   $outputpath= "\./";
}
unlink "$outputpath";
mkdir "$outputpath";
print "the output path is: $outputpath\n";

###################### start to run ######################
open(GENOME, "$options{'genome'}");
system ("cp $options{'genome'} $outputpath/genome.fasta");
my $genomeID; my @totalletter;
while (<GENOME>) {
    if (/>/) {
        s/\>//;
        my @ids=split(/\s+|\t/);
        $genomeID=$ids[0];
    }else{
        push (@totalletter, $_);
    }
}
close GENOME;
my $total=join("",@totalletter);
my @totaleach=split(//,$total);
my $genomelength=grep(/a|t|g|c|u/i,@totaleach);
#print $genomelength,"\n";

open(JBED, "$options{'bed'}");
system ("cp $options{'bed'} $outputpath/junction.bed");
open(GTF, ">$outputpath/genome.gtf");
my $n=0;
while (<JBED>) {
    $n++;
    my @eachbed=split(/\t/);
    print GTF "$genomeID	NCBI	gene	1	$genomelength	.	+	.	gene_id \"subgenome$n\"; gene_name \"$eachbed[0]\"\;\n";
    print GTF "$genomeID	NCBI	exon	1	$eachbed[1]	.	+	.	gene_id \"subgenome$n\"; transcript_id \"Transcript$n\"; gene_name \"$eachbed[0]\"\;\n";
    print GTF "$genomeID	NCBI	exon	$eachbed[2]	$genomelength	.	+	.	gene_id \"subgenome$n\"; transcript_id \"Transcript$n\"; gene_name \"$eachbed[0]\"\;\n";
}
close JBED; close GTF;

system ("hisat2_extract_splice_sites.py $outputpath/genome.gtf > $outputpath/genome_splicesites.txt");
system ("paftools.js gff2bed $outputpath/genome.gtf > $outputpath/genome.bed");

