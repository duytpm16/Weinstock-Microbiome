week=(6 17 24)

for i in {1..19}
do
  for j in "${week[@]}"
  do
    echo "#PBS -l nodes=1:ppn=4,mem=16gb
    module load R/3.5.1

    Rscript single_sample_summary.R $i $j" >> single_summary_chr_${i}_week_${j}.sh
    qsub single_summary_chr_${i}_week_${j}.sh
  done
done
