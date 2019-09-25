#PBS -l nodes=1:ppn=2


for(( i=400; i<=900; i=i+5))
do
  echo " #PBS -l nodes=1:ppn=2
  module load R/3.5.1

  Rscript get_pileups.R $i $(($i + 4))" >> pileup_$i.sh
  qsub pileup_$i.sh
done

