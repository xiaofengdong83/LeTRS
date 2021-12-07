# LeTRS
LeTRS was implemented in Perl programming language, including a main script for identification of leader-TRS junctions and a script for plotting graphs of the results. It accepts fastq files derived from Illumina paired-end and Nanopore cDNA/direct RNA sequencing, and bam files produced by a splicing alignment method with a SARS-CoV-2 genome. By default, LeTRS analyses SARS-CoV-2 by using 10 known leader-TRS junctions and an NCBI reference genome (NC_045512.2), but the user can also provide customized leader-TRS junctions and SARS-CoV-2 or other coronavirus genomes as a reference.<br>

#### 1.Installation:<br>
**(1)Creating an environment in one step**<br>
wget https://github.com/xiaofengdong83/LeTRS/archive/refs/tags/LeTRS_v2.0.1.tar.gz<br>
tar zxvf LeTRS_v2.0.1.tar.gz<br>
cd LeTRS-LeTRS_v2.0.1<br>
conda env create -f my_environment.yml<br>
source activate LeTRS<br>

**(2)Creating an environment step by step**<br>
Third party tools dependencies:<br>
samtools(>=1.11)<br>
hisat2(=2.1.0)<br>
minimap2(=2.17)<br>
[portcullis](https://github.com/maplesond/portcullis)(>=1.1.2)<br>

All these third party tool dependencies should be exported to PATH, so that LeTRS can find them. We suggest installing the portcullis with conda as below:<br>
conda config --add channels defaults<br>
conda config --add channels bioconda<br>
conda config --add channels conda-forge<br>
conda install portcullis=1.2.2<br>

Perl module dependencies:<br>
Getopt::Long<br>
Bio::SeqIO<br>
File::Basename<br>
List::Compare<br>
List::Util<br>
List::Uniq<br>
Statistics::R<br>

R module dependencies (for plotting):<br>
ggplot2<br>
patchwork<br>

#### 2.Usages<br>
Please see the details of each parameter by:<br>

**Examples:**<br>
(1) To analyse ARTIC V3 Nanopore cDNA sequencing data and extract the reads containing the identified leader-TRS junctions in fasta format for pool 1 of amplicon primers (The ARTIC primer_bed can be found in the "primer_bed" folder):<br>
*perl LeTRS.pl -t 16 -extractfasta -pool 1 -Rtch cDNA -mode nanopore -fa example.fastq.gz -primer_bed path_to_primer_V3.bed -o LeTRS_output*<br>

(2) If "-TRSLindependent" option is added, LeTRS will also identify the TRS Leader independent fusion sites in the reads for pool 1 and pool 2 of amplicon primers (The ARTIC primer_bed can be found in the "primer_bed" folder):<br>
*perl LeTRS.pl -t 16 -pool 0 -Rtch cDNA -mode nanopore -TRSLindependent -fa example.fastq.gz -primer_bed path_to_primer_V3.bed -o LeTRS_output*<br>

(3) To analyse direct RNA Nanopore sequencing data and extract the reads containing the identified leader-TRS junctions in fasta format:<br>
*perl LeTRS.pl -t 16 -extractfasta -Rtch RNA -mode nanopore -fq example.fastq.gz -o LeTRS_output*<br>

(4) To analyse direct RNA Nanopore sequencing data with customized leader-TRS junctions and SARS-CoV-2 or other coronavirus genomes as a reference, and extract the reads containing the identified leader-TRS junctions in fasta format (The instruction of making a reference_folder could be found in "readme.txt" of "making_reference_folder_example" folder):<br>
*perl LeTRS.pl -t 16 -extractfasta -Rtch cDNA -mode nanopore -fq example.fastq.gz -primer_bed path_to_custom_primer.bed -o LeTRS_output -ref reference_folder*<br>

(5) To analyse paired end Illumina sequencing data and extract the read pairs containing the identified leader-TRS junctions in fasta format for pool 1 and pool 2 of amplicon primers (The ARTIC primer_bed can be found in the "primer_bed" folder):<br>
*perl LeTRS.pl -t 16 -extractfasta -pool 0 -mode illumina -fq #1.fastq.gz:#2.fastq.gz -primer_bed path_to_primer_V3.bed -o LeTRS_output*<br>

(6) To analyse customized bam file reads derived from any platform aligned by using a splicing mapping method.<br>
*perl LeTRS.pl -t 16 -extractfasta -mode illumina -bam example.bam -o LeTRS_output*<br>

#### 3. Results<br>
The results can be found under the "results" folder in output path, with four tables: known_junction.tab, known_junction_details.tab, novel_junction.tab and novel_junction_details.tab.<br>

**known_junction.tab**<br>
The LeTRS output table for known subgenomic mRNA in the sequencing data. "ref_leader_end" and "peak_leader_end" point to the reference position of the end of the leader and the position of the end of the leader identified in the most common reads (peak count) on the reference genome, and "ref_TRS_start" and "peak_TRS_start" refer to the reference position of the start of the TRS and the position of the start of the TRS identified in the most common reads (peak count) on the reference genome.<br>

**known_junction_details.tab**<br>
The LeTRS output table for details of known subgenomic mRNA in the sequencing data. "peak_leader" and "peak_TRS_start" point to the leader-TRS junctions in known_junction.tab, "ACGAAC" indicates if there is an ACGAAC sequence in the "TRS_seq" (TRS sequences), "20_leader_seq" refers to the 20 nucleotides before the end of the leader, and "ATG_postion" and "first_orf_aa" refer to the first AUG position and translated orf of the sgmRNA.<br>

**novel_junction.tab**<br>
The LeTRS output table for novel subgenomic mRNA in the sequencing data. "leader_end" and "TRS_start" refer to the position of the end of the leader and the position of the start of the TRS identified in the reads >10.<br>

**novel_junction_details.tab**<br>
The LeTRS output table for details of novel subgenomic mRNA in the sequencing data. "peak_leader" and "peak_TRS_start" point to the leader-TRS junctions in novel_junction.tab, "ACGAAC" indicates if there is an ACGAAC sequence in the "TRS_seq" (TRS sequences), "20_leader_seq" refers to the 20 nucleotides before the end of the leader, "AUG_postion" and "first_orf_aa" refer to the first AUG position and translated orf of the sgmRNA, and "known_AUG" indicates if the first AUG position is the same as a known sgmRNA.<br>

**TRS_L_independent_junction.tab**<br>
If "-TRSLindependent" option is added, LeTRS will also identify the TRS Leader independent fusion sites in the reads.<br>

#### 4. Plotting<br>
There is also a perl script that can plot a diagram for the output of LeTRS.pl.<br>

**Examples:**<br>
(1) Plotting the value in the column of "peak_count" in "known_junction.tab" or "nb_count" in the "novel_junction.tab. "-count 1" indicates the first number of each row in the column and "-count 2" indicates the second number of each row in the column, and so on.<br>
*perl LeTRS-plot.pl -count 1 -i known_junction.tab*<br>

(2) Plotting the value in the column of "peak_peak_count_ratio" in "known_junction.tab" or "count_ratio" in the "novel_junction.tab. "-count 1" indicates the first number of each row in the column and "-count 2" indicates the second number of each row in the column, and so on.<br>
*perl LeTRS-plot.pl -ratio 1 -i known_junction.tab*<br>

#### 4. Customized leader-TRS junctions and SARS-CoV-2 or other coronavirus genomes as reference sequences.<br> 
Please the see the "readme.txt" file in the "making_reference_folder_example" folder.<br>
