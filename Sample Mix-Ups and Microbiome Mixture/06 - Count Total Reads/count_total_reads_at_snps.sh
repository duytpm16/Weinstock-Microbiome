week=(6 17 24)


for i in "${week[@]}"
do
  echo "#PBS -l nodes=1:ppn=1
  module load R/3.5.1

  Rscript count_total_reads_at_snps.R $i" >> total_reads_$i.sh
  qsub total_reads_$i.sh
done
