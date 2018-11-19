#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=05:00:00
module load bowtie2/2.3.1







### Directory where .fa genome and bowtie indexes will be stored
mkdir mm10
cd /home/phamd/mm10/
DIR="/home/phamd/mm10/"












### Create sub-directories for chromosomes 1-19 and X in the GRCm38_Index directory
for i in {1..19};
do
  mkdir "${DIR}Chromosome_$i";
done;

mkdir "${DIR}Chromosome_X"









### Moving .fa files to their corresponding directory
cd /home/phamd/
for i in {1..19};
do
  mv "chr${i}.fa" "${DIR}Chromosome_$i";
done;
mv chrX.fa "${DIR}Chromosome_X"









### Create bowtie1 index for the autosome and X chromosome
for i in {1..19};
do
  CHR_DIR="$DIR/Chromosome_$i"
  cd ${CHR_DIR}
  bowtie2-build "chr${i}.fa" "mm10_chr$i";
  cd ${DIR}
done;


CHR_DIR="$DIR/Chromosome_X"
cd ${CHR_DIR}
bowtie2-build "chrX.fa" "mm10_chrX"

