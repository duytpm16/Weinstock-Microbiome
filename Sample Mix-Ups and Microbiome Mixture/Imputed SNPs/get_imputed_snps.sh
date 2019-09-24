


for i in {1..19};
do

  echo "
  #PBS -l nodes=1:ppn=4
  module load R/3.5.1

  Rscript get_imputed_snps.R ${i} 4" >> get_imputed_snps_chr_$i.sh
  
  qsub get_imputed_snps_chr_$i.sh

done 
