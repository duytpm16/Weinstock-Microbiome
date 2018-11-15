#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=05:00:00
module load bowtie/1.1.2







### Directory where .fa genome and bowtie indexes will be stored
cd /home/phamd/GRCm38_Index/
DIR="/home/phamd/GRCm38_Index/"












### Create sub-directories for chromosomes 1-19 and X in the GRCm38_Index directory
for i in {1..19};
do
  mkdir "${DIR}_$i";
done;

mkdir "${DIR}_X"









### Create bowtie1 index for the autosome and X chromosome
for i in {1..19};
do
  CHR_DIR="$DIR/Chromosome_$i"
  cd ${CHR_DIR}
  bowtie-build "Mus_musculus.GRCm38.dna.chromosome.$i.fa" "GRCm38_Chr$i";
  cd ${DIR}
done;


CHR_DIR="$DIR/Chromosome_X"
cd ${CHR_DIR}
bowtie-build "Mus_musculus.GRCm38.dna.chromosome.X.fa" "GRCm38_ChrX"
