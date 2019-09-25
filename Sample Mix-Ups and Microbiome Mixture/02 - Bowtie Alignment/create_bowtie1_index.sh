#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=05:00:00
module load bowtie/1.1.2





### Change to directory containing mm10 fasta files downloaded from:
#     http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
cd /home/phamd/chromFa/






### Concatenate fasta files. Only getting chromosomes 1-19.
fasta_files=""
for i in {1..19};
do
  fasta_files="${fasta_files}chr${i}.fa,";
done








### Seeing if the correct file names were concatenated
echo $fasta_files








### Build mm10 index
bowtie-build $fasta_files mm10.bwt1.autosome
