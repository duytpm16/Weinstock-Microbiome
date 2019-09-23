#PBS -l nodes=1:ppn=2
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1


viewer_data='weinstock_pomp_16s_microbiome_qtl_viewer_v4.RData'
scan1_mat='otu_w6_additive_scan.rds'
dataset_expr='dataset.otu.w6|data|rz'
thr='6'
num_cores='2'
type_scan='additive'
drop='.95'
cis_threshold='NA'
int_mat='NA'



Rscript qtl2_findpeaks.R "viewer_data=$viewer_data" "scan1_mat=$scan1_mat" "dataset_expr=$dataset_expr" "thr=$thr" "num_cores=$num_cores" "type_scan=$type_scan" "drop=$drop" "cis_threshold=$cis_threshold" "int_mat=$int_mat"

