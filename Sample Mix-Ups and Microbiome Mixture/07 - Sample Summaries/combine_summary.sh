#PBS -l nodes=1:ppn=1,walltime=01:00:00

module load R/3.5.1


Rscript combine_summary.R /home/phamd/Weinstock/ sample 6
Rscript combine_summary.R /home/phamd/Weinstock/ sample 17
Rscript combine_summary.R /home/phamd/Weinstock/ sample 24



Rscript combine_summary.R /home/phamd/Weinstock/ paired 6
Rscript combine_summary.R /home/phamd/Weinstock/ paired 17
Rscript combine_summary.R /home/phamd/Weinstock/ paired 24
