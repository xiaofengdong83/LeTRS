# LeTRS
LeTRS was implemented in Perl programming language, including a main script for identification of leader-TRS junctions and a script for plotting graphs of the results. It accepts fastq files derived from illumina paired-end and nanopore cDNA/direct RNA sequencing, and bam files produced by a splicing alignment method with a sars-cov-2 genome. By default, LeTRS analyses the sars-cov-2 by using 10 known leader-TRS junctions and an NCBI reference genome (NC_045512.2), but user can also provide customize leader-TRS junctions and sars-cov-2 or other coronavirus genome as reference.<br>

#### 1.Installation:<br>
**Three party tools dependences:**<br>
samtools (>=1.11)<br>
hisat2(>=2.1.0)<br>
minimap2(>=2.17)<br>
[portcullis](https://github.com/maplesond/portcullis)(=1.1.2)<br>

All these three party tools dependences should be exported to PATH, so that LsTRS can din them. We suggest to install the portcullis with conda as below:<br>
conda config --add channels defaults<br>
conda config --add channels bioconda<br>
conda config --add channels conda-forge<br>
conda install portcullis=1.1.2<br>

**Perl module dependences:**<br>
Getopt::Long<br>
Bio::SeqIO<br>
File::Basename<br>
List::Compare<br>
List::Util<br>
List::Uniq<br>
Statistics::R<br>

**R module dependences (for plotting):**<br>
ggplot2<br>

#### 2.Usages<br>
Please see the details of each parameters by:<br>

**Examples:**<br>
(1) To analyse ARTIC V3 noropore cDNA sequencing data and extract the reads contain the identified leader-TRS junctions in fasta format (The ARTIC primer_bed can be found in the "primer_bed" folder):<br>
*perl LeTRS.pl -t 16 -extractfasta -Rtch cDNA -mode noropore -fa example.fastq.gz -primer_bed primer_V3.bed -o LeTRS_output*<br>

(2) To analyse direct RNA noropore sequencing data and extract the reads contain the identified leader-TRS junctions in fasta format (The ARTIC primer_bed can be found in the "primer_bed" folder):<br>
*perl LeTRS.pl -t 16 -extractfasta -Rtch RNA -mode noropore -fq example.fastq.gz -primer_bed primer_V3.bed -o LeTRS_output*<br>

(3) To analyse direct RNA noropore sequencing data with customize leader-TRS junctions and sars-cov-2 or other coronavirus genome as reference, and extract the reads contain the identified leader-TRS junctions in fasta format (The ARTIC primer_bed can be found in the "primer_bed" folder):<br>
*perl LeTRS.pl -t 16 -extractfasta -Rtch RNA -mode noropore -fq example.fastq.gz -primer_bed primer_V3.bed -o LeTRS_output -ref reference_folder*<br>

(4) To analyse paired end illumina sequencing data and extract the reads contain the identified leader-TRS junctions in fasta format (The ARTIC primer_bed can be found in the "primer_bed" folder):<br>
*perl LeTRS.pl -t 16 -extractfasta -mode illumia -fq #1.fasq.gz:#2.fasq.gz -primer_bed primer_V3.bed -o LeTRS_output*<br>

(5) To analyse customize bam file reads derived from any platform aligned by using a splicing mapping method.<br>
*perl LeTRS.pl -t 16 -extractfasta -mode illumia -bam example.bam -o LeTRS_output*<br>

#### 3. Results<br>
The results can be found under the "results" folder in output path, with for tables: known_junction.tab, known_junction_details.tab, novel_junction.tab and novel_junction_details.tab.<br>

**known_junction.tab**<br>
The LeTRS output table for known subgenomic mRNA in the sequencing data. "ref_leader_end" and "peak_leader_end" point to the reference position of the end of leader and the position of the end of leader identified in the most common reads (peak count) on the reference genome, and "ref_TRS_start" and "peak_TRS_start" refer to the reference position of the start of TRS and the position of the start of TRS identified in the most common reads (peak count) on the reference genome.<br>

**known_junction_details.tab**<br>
The LeTRS output table for details of known subgenomic mRNA in the sequencing data. "peak_leader" and "peak_TRS_start" point to the leader-TRS junctions in known_junction.tab, "ACGAAC" indicates if there is a ACGAAC sequence in the "TRS_seq" (TRS sequences), "20_leader_seq" refers to the 20 nucleotides before the end of leader, and "ATG_postion" and "first_orf_aa" refer to the first AUG position and translated orf of the sgmRNA.<br>

**novel_junction.tab**<br>
The LeTRS output table for novel subgenomic mRNA in the sequencing data. "leader_end" and "TRS_start" refer to the position of the end of leader and the position of the start of TRS identified in the reads >10.<br>

**novel_junction_details.tab**<br>
The LeTRS output table for details of novel subgenomic mRNA in the sequencing data. "peak_leader" and "peak_TRS_start" point to the leader-TRS junctions in novel_junction.tab, "ACGAAC" indicates if there is a ACGAAC sequences in the "TRS_seq" (TRS sequences), "20_leader_seq" refers to the 20 nucleotides before the end of leader, "AUG_postion" and "first_orf_aa" refer to the first AUG position and translated orf of the sgmRNA, and "known_AUG" indicates if the first AUG position same as a known sgmRNA.<br>

#### 4. Plotting<br>
There is also a perl script that can plot a diagram for the output of LeTRS.pl.<br>

**Examples:**<br>
(1) Plotting the value in the column of "peak_count" in "known_junction.tab" or "nb_count" in the "novel_junction.tab. "-count 1" indicates the first number of each row in the column and "-count 2" indicates the second number of each row in the column, and so on.<br>
*perl LeTRS-plot.pl -count 1 -i known_junction.tab*<br>

(2) Plotting the value in the column of "peak_peak_count_ratio" in "known_junction.tab" or "count_ratio" in the "novel_junction.tab. "-count 1" indicates the first number of each row in the column and "-count 2" indicates the second number of each row in the column, and so on.<br>
*perl LeTRS-plot.pl -ratio 1 -i known_junction.tab*<br>

#### 4. Customize leader-TRS junctions and sars-cov-2 or other coronavirus genome as reference.<br>
Please the see the "readme.txt" file in the "making_reference_folder_example" folder.<br>
