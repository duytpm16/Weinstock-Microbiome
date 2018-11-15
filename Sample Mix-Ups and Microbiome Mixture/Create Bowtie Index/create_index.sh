#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=05:00:00
module load bowtie/1.1.2




cd /home/phamd/GRCm38_Index/
DIR="/home/phamd/GRCm38_Index/"




### Create bowtie1 index for the autosome
for i in {1..19};
do
  CHR_DIR="$DIR/Chromosome_$i"
  cd ${CHR_DIR}
  bowtie-build "Mus_musculus.GRCm38.dna.chromosome.$i.fa" "GRCm38_Chr$i";
  cd ${DIR}
done;






### Create bowtie index for the X chromosome
CHR_DIR="$DIR/Chromosome_X"
cd ${CHR_DIR}
bowtie-build "Mus_musculus.GRCm38.dna.chromosome.X.fa" "GRCm38_ChrX"
