### Shell script to run qtl2 on HPC with data in QTL viewer format.
#   This is for additive scan (int_name = 'NA').
#   I did not run this in chunks, hence the 'NA'.
#   See https://github.com/duytpm16/qtl2-HPC-pipeline for input details.

#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00

module load R/3.5.1

viewer_data='weinstock_pomp_16s_microbiome_qtl_viewer_v4.RData'
dataset_expr='dataset.otu.w24|data|rz'
num_cores='8'
int_name='NA'
chunk_number='NA'
chunk_size='NA'



Rscript qtl2_scan1.R viewer_data=$viewer_data dataset_expr=$dataset_expr num_cores=$num_cores int_name=$int_name chunk_number=$chunk_number chunk_size=$chunk_size
