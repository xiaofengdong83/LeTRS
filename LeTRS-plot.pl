#!/usr/bin/perl -w
use strict;  
use Getopt::Long;
use List::Uniq ':all';
use Statistics::R;

my %options;
my @standard_options =("help|h!",
                       "i=s",
                       "o=s",
                       "count=s",
                       "ratio=s"
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
  perl LeTRS-plot.pl -count 1 -i known_junction.tab
  perl LeTRS-plot.pl -ratio 1 -i known_junction.tab
   
  -i                input files \"known_junction.tab\" or \"novel_junction.tab\".
  -o                output path, \"./\" by default.
  -count            1,2,3.. indicates the value in the column of peak_count in \"known_junction.tab\"
                    or nb_count in the \"novel_junction.tab\".
  -ratio            1,2,3.. indicates the value in the column of peak_count_ratio in \"known_junction.tab\"
                    or count_ratio in the \"novel_junction.tab\".
                    
  -h/-help          Produce help message.\n\n";
    exit;
}

my $checkcount; my $checkratio;
if (exists $options{'count'}) {
   $checkcount=1;
}elsif(!exists $options{'count'}) {
   $checkcount=0;
}
if (exists $options{'ratio'}) {
   $checkratio=1;
}elsif(!exists $options{'ratio'}) {
   $checkratio=0;
}
if ($checkcount==1 and $checkratio==1) {
    print "please only provide one of \"-count\" and \"-ratio\"\n";
    exit;
}
if ($checkcount==0 and $checkratio==0) {
    print "please only provide \"-count\" or \"-ratio\"\n";
    exit;
}
if ($checkcount==1 and $checkratio==0) {
    print "the value is: count $options{'count'}\n";
}
if ($checkcount==0 and $checkratio==1) {
    print "the value is: ratio $options{'ratio'}\n";
}

if (!exists $options{'i'}) {
   print "please provide input table\n";
   exit;
}
print "the input path is: $options{'i'}\n";

my $outputpath;
if (exists $options{'o'}) {
   $outputpath=$options{'o'};  
}else{
   $outputpath= "\./";
}
print "the output path is: $outputpath\n";



###################### start to run ######################
open(TABLE,"$options{'i'}");
open(TABLERR,">$outputpath/known_junction_tmp2.tab");
print TABLERR "leaderorf\tcount\tlevel\n";
my @tables=<TABLE>;
my @collection; my @leaders; my @leaders2; my @orfs; my $numfortmp1;
foreach (@tables) {
   if ($tables[0]=~/^subgenome\tref_leader_end/) {
      unless (/^subgenome\tref_leader_end|The numbers in the bracket|Normalized count\=|Total number of read mapped/) {
         my @each1=split(/\t/);
         my @each2;
         if (exists $options{'count'}) {
            @each2=split(/\(|\,/,$each1[5]);
         }
         if (exists $options{'ratio'}) {
            @each2=split(/\(|\,/,$each1[6]);
         }
         
         unless($each1[2]==0) {
            my $plotvalue;
            if (exists $options{'count'}) {
                $plotvalue=$options{'count'}-1;
            }
            if (exists $options{'ratio'}) {
                $plotvalue=$options{'ratio'}-1;
            }
            #print "$each1[0]\t$each1[1]\t$each1[3]\t$each1[5]\t$each2[$plotvalue]\n";
            push (@collection, "leader\_$each1[0]\_$each1[2]\tleader\_$each1[2]\t$each1[2]\n");
            push (@collection, "leader\_$each1[0]\_$each1[2]\t$each1[0]\t$each1[0]\n");
            $numfortmp1++;
            print TABLERR "leader\_$each1[0]\_$each1[2]\t$each2[$plotvalue]\t$numfortmp1\n";
            push (@leaders, $each1[2]);
            push (@leaders2, "leader\_$each1[0]\_$each1[2]");
            push (@orfs, $each1[0]);
         }
      }
   }
   
   if ($tables[0]=~/^subgenome\tleader_end/) {
      unless (/^subgenome\tleader_end|The numbers in the bracket|Normalized count\=|Total number of read mapped/) {
         my @each1=split(/\t/);
         my @each2;
         if (exists $options{'count'}) {
            @each2=split(/\(|\,/,$each1[3]);
         }
         if (exists $options{'ratio'}) {
            @each2=split(/\(|\,/,$each1[4]);
         }
         
         my $plotvalue;
         if (exists $options{'count'}) {
             $plotvalue=$options{'count'}-1;
         }
         if (exists $options{'ratio'}) {
             $plotvalue=$options{'ratio'}-1;
         }
            
         #print "$each1[0]\t$each1[1]\t$each1[2]\t$each2[$plotvalue]\n";
         push (@collection, "$each1[1]\_$each1[2]\tleader\_$each1[1]\t$each1[1]\n");
         push (@collection, "$each1[1]\_$each1[2]\t$each1[2]\t$each1[2]\n");
         $numfortmp1++;
         print TABLERR "$each1[1]\_$each1[2]\t$each2[$plotvalue]\t$numfortmp1\n";
         push (@leaders, "$each1[1]");
         push (@leaders2, "$each1[1]\_$each1[2]");
         push (@orfs, $each1[2]);
      }
   }
   
}
close TABLE;

my %hashleader; my $numleaders=0;
my @uniqleaders=uniq(@leaders);
foreach (@uniqleaders) {
   $numleaders++;
   $hashleader{$_}=$numleaders;
}

my %hashleader2; my $numleaders2=0;
foreach (@leaders2) {
   $numleaders2++;
   $hashleader2{$_}=$numleaders2;
}

my %hashorf; my $numorf=$numleaders;
foreach (@orfs) {
   $numorf++;
   $hashorf{$_}=$numorf;
}

open(TABLER,">$outputpath/known_junction_tmp1.tab");
print TABLER "leader\tjunction\tlabel\tlevels\tlevels2\n";
foreach (@collection) {
   chomp;
   my @eachcollection=split(/\t/);
   
   if ($eachcollection[1]=~/leader_/) {
      print TABLER "$eachcollection[0]\t$eachcollection[1]\t$eachcollection[2]\t$hashleader{$eachcollection[2]}\t$hashleader2{$eachcollection[0]}\n";
   }else{
      print TABLER "$eachcollection[0]\tTRS_$eachcollection[1]\t$eachcollection[2]\t$hashorf{$eachcollection[2]}\t$hashleader2{$eachcollection[0]}\n";
   }
}
close TABLER;

my $R = Statistics::R->new();
my $input_value1 = "$outputpath/known_junction_tmp1.tab";
my $input_value2 = "$outputpath/known_junction_tmp2.tab";
my $input_value3 = "$outputpath/leader-TRS.pdf";
my $input_value4 = $#orfs+1;
$R->set('tmp1', $input_value1);
$R->set('tmp2', $input_value2);
$R->set('outputplot', $input_value3);
$R->set('plotsize', $input_value4);
if (exists $options{'count'}) {
   $R->run('library(ggplot2)
   library(patchwork)

   plotdatabarchart<-read.table(tmp2, head=T, row.names = NULL)
   plotdatabarchart$leaderorf <- factor(plotdatabarchart$leaderorf,levels= unique(plotdatabarchart[order(plotdatabarchart$level), "leaderorf"]))

   # Basic barplot
   plot1<-ggplot(data= plotdatabarchart, aes(x=leaderorf, y=count)) + geom_bar(stat="identity", fill="red", width=0.3)+ geom_text(aes(y= count, label= count), vjust=-0.1, color="black")+ labs(y = "Count")+ theme(panel.background = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank())

   plotdata<-read.table(tmp1, head=T, row.names = NULL)

   plotdata$junction <- factor(plotdata$junction,levels= rev(unique(plotdata[order(plotdata$levels), "junction"])))
   plotdata$leader <- factor(plotdata$leader,levels= unique(plotdata[order(plotdata$levels2), "leader"]))

   plot2<-ggplot(data= plotdata, aes(x= leader,y= junction)) + geom_line(aes(group = leader))+ geom_point(color= "blue", size=2) + labs(y = "Position on the reference genome (nt)")+ theme(axis.line = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank())

   #plotsizewith<-plotsize
   #plotsizeheight<-plotsize/1.5
   #pdf(file = outputplot, width = plotsizewith, height = plotsizeheight)
   pdf(file = outputplot)
   plot1 + plot2 + plot_layout(ncol = 1, heights = c(6, 6))
   dev.off()');
}

if (exists $options{'ratio'}) {
$R->run('library(ggplot2)
   library(patchwork)

   plotdatabarchart<-read.table(tmp2, head=T, row.names = NULL)
   plotdatabarchart$leaderorf <- factor(plotdatabarchart$leaderorf,levels= unique(plotdatabarchart[order(plotdatabarchart$level), "leaderorf"]))

   # Basic barplot
   plot1<-ggplot(data= plotdatabarchart, aes(x=leaderorf, y=count)) + geom_bar(stat="identity", fill="red", width=0.3)+ geom_text(aes(y= count, label= count), vjust=-0.1, color="black")+ labs(y = "Normalized count")+ theme(panel.background = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank())

   plotdata<-read.table(tmp1, head=T, row.names = NULL)

   plotdata$junction <- factor(plotdata$junction,levels= rev(unique(plotdata[order(plotdata$levels), "junction"])))
   plotdata$leader <- factor(plotdata$leader,levels= unique(plotdata[order(plotdata$levels2), "leader"]))

   plot2<-ggplot(data= plotdata, aes(x= leader,y= junction)) + geom_line(aes(group = leader))+ geom_point(color= "blue", size=2) + labs(y = "Position on the reference genome (nt)")+ theme(axis.line = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank())
   
   #plotsizewith<-plotsize
   #plotsizeheight<-plotsize/1.5
   #pdf(file = outputplot, width = plotsizewith, height = plotsizeheight)
   pdf(file = outputplot)
   plot1 + plot2 + plot_layout(ncol = 1, heights = c(6, 6))
   dev.off()');
}
unlink "$outputpath/known_junction_tmp1.tab";
unlink "$outputpath/known_junction_tmp2.tab";

