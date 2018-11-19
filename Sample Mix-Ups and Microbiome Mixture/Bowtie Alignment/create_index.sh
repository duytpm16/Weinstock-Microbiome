#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=05:00:00
module load bowtie2/2.3.1


cd chromFA/

# Concatenate fasta files
fasta_files=""
for i in {1..19};
do
  fasta_files="${fasta_files},chr${i}.fa";
done

fasta_files="{$fasta_files},chrX.fa"





echo $fasta_files




bowtie2-build fasta_files
