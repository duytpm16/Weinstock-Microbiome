### Options
options(stringsAsFactors = FALSE)





### Command line arguments / variables to change
args = commandArgs(trailingOnly = TRUE)
chr  <- args[1]
week <- args[2]
snps_directory  <- "/projects/churchill-lab/data/Weinstock/Pomp_Benson/genotypes/imputed_snps/"
fastq_directory <- "/projects/churchill-lab/data/Weinstock/Pomp_Benson/fastq/"
save_directory  <- "/home/phamd/Weinstock/"









### Load the corresponding imputed SNPs and SNP info for that chromosome
load(paste0(snps_directory, "imp_snp_", chr, ".Rdata"))
snpinfo$pos_bp <- round(snpinfo$pos * 1e6) 
save_directory <- paste0(save_directory,'week_',week,'/') 










### Change directory to where fastq are store and get sample ids
setwd(fastq_directory)
sample_num   <- sapply(dir()[grep('MGS_DO', dir())], function(s) unlist(strsplit(s, '_'))[3])
sample_exist <- sapply(dir()[grep('MGS_DO', dir())], function(s) dir.exists(paste0(s,'/Week_',week,'/')))
stopifnot(names(sample_num) == names(sample_exist))
sample_num   <- sample_num[sample_exist]

sample_id <- paste0('DP.DO2.',sample_num,'.F')







### Creating object to store results
sample_results <- vector("list", length(sample_num))
names(sample_results) <- sample_id












### Loop through each MGS_DO_* directory
for(index in 1:length(sample_num)){
  
    # Get directory name
    sample <- paste0('MGS_DO_', sample_num[index])
    sample_week_directory <- paste0(fastq_directory, sample, '/Week_', week,'/bowtie1_run')
  
  
  
  
  
    # Change to MGS_DO_* directory and week directory if it exists
    if(dir.exists(sample_week_directory)){
       setwd(sample_week_directory)
       print(sample_week_directory)
       print(sample)  
    
    
    
    
    
    
       # Read in read counts of one chromosome of one sample 
       sample_read_counts <- readRDS(paste0(sample,'_Week_',week,'_readcounts_chr_',chr,'.rds'))
 
    
  
    
    
    
       # Making sure all positions and snp id are the same
       filtered_snpinfo <- snpinfo[snpinfo$snp_id %in% colnames(imp_snps),]
       aligning_pos <- match(sample_read_counts$pos, snpinfo$pos_bp)
       stopifnot(!any(is.na(aligning_pos)))


       imp_snps_col <- filtered_snpinfo$snp_id[aligning_pos]
       print(length(imp_snps_col))

       stopifnot(length(imp_snps_col) == nrow(sample_read_counts))
       stopifnot(all(filtered_snpinfo$pos_bp[aligning_pos] == sample_read_counts$pos))


    
    
    
    
    
       # Create object to contain the results for the single samples
       sample_results[[index]]           <- array(0, dim=c(nrow(imp_snps), 3, 2))
       dimnames(sample_results[[index]]) <- list(rownames(imp_snps), c("AA", "AB", "BB"), c("A", "B"))
    
    
    
       for(i in 1:nrow(imp_snps)) {
           g <- imp_snps[i, imp_snps_col]
        
           for(j in 1:3) {
             
               sample_results[[index]][i,j,1] <- sum(sample_read_counts$n1[!is.na(g) & g==j])
               sample_results[[index]][i,j,2] <- sum(sample_read_counts$n2[!is.na(g) & g==j])
             
           } # for j
       } # for i
    
    
    
    
    
    
  } # If
} # Loop through each MB sample






# Save                         
saveRDS(sample_results, file = paste0(save_directory, "sample_results_week_",week,"_chr_", chr, ".rds"))
