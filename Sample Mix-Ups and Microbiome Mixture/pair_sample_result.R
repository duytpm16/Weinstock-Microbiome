# get summaries for each MB sample against all DNA samples, one chr
args = commandArgs(trailingOnly = TRUE)
home_directory  <- '/home/phamd/'
fastq_directory <- "/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/"
chr             <- args[1]
week            <- args[2]



# load the corresponding imputed SNPs for that chromosome
load(paste0(home_directory,"imputed_snps_chr_", chr, ".RData"))
imp_snps             <- get(paste0('imp_snps_chr', chr))
index_snpinfo        <- get(paste0('indexed_snpinfo_chr', chr))
index_snpinfo$pos_bp <- round(index_snpinfo$pos * 1e6)
save_directory       <- paste0(home_directory,'week_',week,'/')





setwd(fastq_directory)
sample_num     <- sapply(dir(), function(s) unlist(strsplit(s, '_'))[3])
sample_id      <- paste0('DPDP.DO2.',sample_num,'.F')

print(sample_id)




pair_results   <- vector("list", length(dir()))
names(pair_results)   <- sample_id







### Loop through each MGS_DO_* directory
for(index in 1:length(sample_num)){
  
    # Change to MGS_DO_* directory and week directory
    sample <- paste0('MGS_DO_', sample_num[index])
    sample_week_directory <- paste0(fastq_directory, sample, '/Week_', week)
    if(dir.exists(sample_week_directory)){
       setwd(sample_week_directory)
       print(sample_week_directory)
  
  
  
  
       # Read in read counts of one chromosome for one sample 
       sample_read_counts <- readRDS(paste0(sample,'_Week_',week,'_readcounts_chr_',chr,'.rds'))
       sample_read_counts <- sample_read_counts[!(sample_read_counts$Major_Count == 0 & sample_read_counts$Minor_Count == 0),]

  
  
  
    
       aligning_pos           <- sample_read_counts$pos[sample_read_counts$pos %in% index_snpinfo$pos_bp]
       filtered_index_snpinfo <- index_snpinfo[index_snpinfo$pos_bp %in% aligning_pos,]
       sample_read_counts     <- sample_read_counts[sample_read_counts$pos %in% aligning_pos,]
       imp_snps_col           <- colnames(imp_snps)[colnames(imp_snps) %in% filtered_index_snpinfo$snp]
  
  
       print(length(imp_snps_col))    
    
  
       # create object to contain the results for sample pairs
       pair_results[[index]] <- array(0, dim=c(nrow(imp_snps), 3, 3, 2))
       dimnames(pair_results[[index]]) <- list(rownames(imp_snps), c("AA", "AB", "BB"),
                                            c("AA", "AB", "BB"), c("A", "B"))
  
       g0 <- imp_snps[sample_id[index], imp_snps_col]
       for(i in 1:nrow(imp_snps)) {
    
           g <- imp_snps[i, imp_snps_col]
     
           for(j in 1:3) {
      
               for(k in 1:3){
                  
                   pair_results[[index]][i, j, k, 1] <- sum(sample_read_counts$Major_Count[!is.na(g0) & g0==j & !is.na(g) & g==k])
                   pair_results[[index]][i, j, k, 2] <- sum(sample_read_counts$Minor_Count[!is.na(g0) & g0==j & !is.na(g) & g==k])
              
               } # for k
              
           } # for j
    
       } # for i 
  
    } # if  
    
} # loop over MB samples



saveRDS(pair_results, paste0(save_directory, "paired_results_week_",week,"_chr_", chr, ".rds"))
