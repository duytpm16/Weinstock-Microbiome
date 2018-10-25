### HPC shell script for qsubing qtl2_scan1_pomp_microbiome.R script
#    Ran 1 at a time by commenting out the others and qsubing

### Parameters after Rscript:
#       1.) Name of R script file
#       2.) Which taxa
#       3.) Number of cores to run
#       4.) Logical value: TRUE = qtl scan with rankZ transform data
#                          FALSE = qtl scan with normalized data

#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00


module load R/3.5.1

Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome otu 8 TRUE
Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome genus 8 TRUE
Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome family 8 TRUE
Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome class 8 TRUE
Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome order 8 TRUE
Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome phylum 8 TRUE

