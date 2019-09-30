#PBS -l nodes=1:ppn=8

module load R/3.5.1


Rscript pair_mbmixture_analysis.R /home/phamd/Weinstock/paired_results_week_6_all.rds 6 8 
Rscript pair_mbmixture_analysis.R /home/phamd/Weinstock/paired_results_week_17_all.rds 17 8
Rscript pair_mbmixture_analysis.R /home/phamd/Weinstock/paired_results_week_24_all.rds 24 8
