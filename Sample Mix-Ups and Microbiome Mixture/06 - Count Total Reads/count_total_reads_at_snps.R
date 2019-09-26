### Options
options(stringsAsFactors = FALSE)





### Command line arguments / variables to change
args <- commandArgs(trailingOnly = TRUE)
week <- args[1]
fastq_directory <- "/projects/churchill-lab/data/Weinstock/Pomp_Benson/fastq/"
save_directory  <- "/home/phamd/Weinstock/"



print(week)







### Change directory to where fastq are store and get sample ids
setwd(fastq_directory)
sample_num   <- sapply(dir()[grep('MGS_DO', dir())], function(s) unlist(strsplit(s, '_'))[3])
sample_exist <- sapply(dir()[grep('MGS_DO', dir())], function(s) dir.exists(paste0(s,'/Week_',week,'/')))
stopifnot(names(sample_num) == names(sample_exist))
sample_num   <- sample_num[sample_exist]

sample_id <- paste0('DP.DO2.',sample_num,'.F')







### Creating object to store results
readcounts <- matrix(nrow=length(sample_num), ncol=19)
dimnames(readcounts) <- list(sample_id, 1:19)







### Loop through each MGS_DO_* directory
for(index in 1:length(sample_num)){
  
    # Get directory name
    sample <- paste0('MGS_DO_', sample_num[index])
    sample_week_directory <- paste0(fastq_directory, sample, '/Week_', week,'/bowtie1_run_v2/')
    setwd(sample_week_directory)

    for(chr in 1:19){
        rd <- readRDS(grep(paste0('readcounts_chr_',chr,'.rds'), dir(), value = TRUE))
        readcounts[sample_id[index], chr] <- sum(rd$n1 + rd$n2)
    }

}




saveRDS(readcounts, file = paste0(save_directory,'total_readcounts_week_', week, '.rds'))
  
  
  
