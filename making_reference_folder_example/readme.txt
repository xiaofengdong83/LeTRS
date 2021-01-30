making_ref.pl can help user to create reference folder for customize leader-TRS junctions and sars-cov-2 or other coronavirus genome as reference.

The input files are fasta file of a coronavirus reference genome and a four columns bed file including subgenomic name, position of the end of leader, position of the start of TRS and position of the end of TRS.

to run this below as a example:
perl making_ref.pl -genome genome.fasta -bed junction.bed -o reference_folder
