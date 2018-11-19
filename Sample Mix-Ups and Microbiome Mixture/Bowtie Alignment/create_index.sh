#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=05:00:00
module load bowtie2/2.3.1



cd chromFA/

# Concatenate fasta files
fasta_files=""
for i in {1..19};
do
  fasta_files="${fasta_files}chr${i}.fa,";
done

fasta_files="${fasta_files}chrX.fa"







# Seeing if the correct file names were concatenated
echo $fasta_files








# Build mm10 index
bowtie2-build $fasta_files mm10.genome
