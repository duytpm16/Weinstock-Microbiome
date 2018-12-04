#!/bin/bash
#PBS -q batch


### Parallel imputed snps for each chromosome
#module load R/3.5.1
#for i in {18..19}
#do
#echo "#PBS -q batch
##PBS -l nodes=1:ppn=8

#module load R/3.5.1
#Rscript get_imputed_snps.R $i" >> get_imputed_snps_$i.sh
#qsub get_imputed_snps_$i.sh
#done

#Rscript get_imputed_snps.R X







### Parrallel get pileups for each sample
#for((n=400; n<=855; n = n+5))
#do

#start=$(($n))
#end=$(($n + 4))


#echo "#PBS -q batch
##PBS -l nodes=1:ppn=4
##PBS -l walltime=24:00:00

#module load R/3.5.1

#Rscript get_pileup.R $start $end" >> get_pileup_$n.sh
#qsub get_pileup_$n.sh

#done







## Parrallel get readcounts for each sample
#for((n=700; n<=705; n = n+5))
#do

#start=$(($n))
#end=$(($n + 4))


#echo "#PBS -q batch
##PBS -l nodes=1:ppn=4
##PBS -l walltime=24:00:00

#module load R/3.5.1
#
#Rscript get_readcounts.R $start $end" >> get_readcounts_$n.sh
#qsub get_readcounts_$n.sh

#done







### Parallel imputed snps for each chromosome
week=("6" "17" "24")
chromosome=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X")
#for w in "${week[@]}"
#do
#  for i in "${chromosome[@]}"
#  do
#    echo "#PBS -q batch
#    #PBS -l nodes=1:ppn=4

#    module load R/3.5.1
#    Rscript single_sample_result.R $i $w" >> "sample_results_wk${w}_${i}.sh"
#    qsub "sample_results_wk${w}_${i}.sh"
#  done
#done






for w in "${week[@]}"
do
  for i in "${chromosome[@]}"
  do
    echo "#PBS -q batch
    #PBS -l nodes=1:ppn=1
    #PBS -l walltime=72:00:00 
   
    module load R/3.5.1
    Rscript pair_sample_result.R $i $w" >> "pair_results_wk${w}_${i}.sh"
    qsub "pair_results_wk${w}_${i}.sh"
  done
done

