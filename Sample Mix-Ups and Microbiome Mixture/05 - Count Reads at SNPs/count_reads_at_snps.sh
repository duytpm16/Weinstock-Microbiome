#PBS -l nodes=1:ppn=1


for(( i=400; i<=900; i=i+5))
do
  echo " #PBS -l nodes=1:ppn=1
  module load R/3.5.1
  Rscript count_reads_at_snps.R $i $(($i + 4))" >> count_reads_$i.sh
  qsub count_reads_$i.sh
done
