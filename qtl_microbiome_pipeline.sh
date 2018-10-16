#PBS -l nodes=1:ppn=32
#PBS -q batch
#PBS -l walltime=72:00:00


module load R/3.5.1


#Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome genus 32 TRUE
#Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome family 32 TRUE
#Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome class 32 TRUE
#Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome order 32 TRUE
#Rscript qtl2_scan1_pomp_microbiome.R weinstock_16s_microbiome phylum 32 TRUE

