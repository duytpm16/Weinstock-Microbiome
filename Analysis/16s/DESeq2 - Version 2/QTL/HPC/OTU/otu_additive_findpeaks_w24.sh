### Shell script to run qtl2 'findpeaks' on HPC with data in QTL viewer format.
#   This is for additive scan (type_scan='additive', int_mat = 'NA').
#   I saved all QTLs with LOD > 6 (thr=6).
#   I used a 95% Bayesian confidence interval to estimate interval of QTLs.
#   See https://github.com/duytpm16/qtl2-HPC-pipeline for input details.

#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1


viewer_data='weinstock_pomp_16s_microbiome_qtl_viewer_v4.RData'
scan1_mat='otu_w24_additive_scan.rds'
dataset_expr='dataset.otu.w24|data|rz'
thr='6'
num_cores='8'
type_scan='additive'
drop='.95'
cis_threshold='NA'
int_mat='NA'



Rscript qtl2_findpeaks.R "viewer_data=$viewer_data" "scan1_mat=$scan1_mat" "dataset_expr=$dataset_expr" "thr=$thr" "num_cores=$num_cores" "type_scan=$type_scan" "drop=$drop" "cis_threshold=$cis_threshold" "int_mat=$int_mat"

