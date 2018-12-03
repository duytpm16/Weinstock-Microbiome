
# get summaries for each MB sample against all DNA samples, one chr
args = commandArgs(trailingOnly = TRUE)
home_directory  <- '/home/phamd/'
fastq_directory <- "/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/"
chr             <- args[1]
week            <- args[2]


# load the corresponding imputed SNPs for that chromosome
load(paste0("/home/phamd/imputed_snps_chr_", chr, ".RData"))
imp_snps             <- get(paste0('imp_snps_chr', chr))
index_snpinfo        <- get(paste0('indexed_snpinfo_chr', chr))
index_snpinfo$pos_bp <- round(index_snpinfo$pos * 1e6)
save_directory       <- paste0(home_directory,'week_',week,'/')





setwd(fastq_directory)
sample_num     <- sapply(dir(), function(s) unlist(strsplit(s, '_'))[3])
sample_id      <- paste0('DPDP.DO2.',sample_num,'.F')






sample_results   <- vector("list", length(dir()))
names(sample_results)   <- sample_id







### Loop through each MGS_DO_* directory
for(index in 1:length(sample_num)){
  
    # Change to MGS_DO_* directory and week directory
    sample <- paste0('MGS_DO_', sample_num[index])
    sample_week_directory <- paste0(fastq_directory, sample, '/Week_', week)
    setwd(sample_week_directory)
    print(sample_week_directory)
    print(sample)  
  
  
  
    # Read in read counts of one chromosome for one sample 
    sample_read_counts <- readRDS(paste0(sample,'_Week_',week,'_readcounts_chr_',chr,'.rds'))
    sample_read_counts <- sample_read_counts[!(sample_read_counts$Major_Count == 0 & sample_read_counts$Minor_Count == 0),] 
  
  
  
  
  
    # Making sure all positions are found
    snpinfo_row <- match(sample_read_counts$pos, index_snpinfo$pos_bp)
    snpinfo_row <- snpinfo_row[complete.cases(snpinfo_row)]
    stopifnot(!any(is.na(sample_read_counts$pos)))
  
  
  
  
    # Extract snp names that match positions found. Will be used to extract the imputed snps matrix columns
    imp_snps_col <- index_snpinfo$snp[snpinfo_row]
    print(length(imp_snps_col)) 
  
  
  
  
    # Create object to contain the results for the single samples
    sample_results[[index]] <- array(0, dim=c(nrow(imp_snps), 3, 2))
    dimnames(sample_results[[index]]) <- list(rownames(imp_snps), c("AA", "AB", "BB"), c("A", "B"))
  
  
  
  
  
    for(i in 1:nrow(imp_snps)) {
        g <- imp_snps[i, colnames(imp_snps) %in% imp_snps_col]
        for(j in 1:3) {
            sample_results[[index]][i,j,1] <- sum(sample_read_counts$Major_Count[!is.na(g) & g==j])
            sample_results[[index]][i,j,2] <- sum(sample_read_counts$Minor_Count[!is.na(g) & g==j])
        }
    }

}
  
  
saveRDS(sample_results, file = paste0(save_directory, "sample_results_week_",week,"_chr_", chr, ".rds"))
