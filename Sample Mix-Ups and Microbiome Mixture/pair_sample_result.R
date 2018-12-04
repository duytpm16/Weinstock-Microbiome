### Options
options(stringsAsFactors = FALSE)







### Command line arguments / variables to change
#     1.) home_directory : directory where my imputed snps .RData file are stored + where I want to save output of this script
#     2.) fastq_directory: directory where each sub-directory is a DO sample and contains the fastq, pileup, and read counts
#     3.) chr            : chromosome readcounts to get
#     4.) week           : which week to get readcounts from
args = commandArgs(trailingOnly = TRUE)
home_directory  <- '/home/phamd/'
fastq_directory <- "/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/"
chr             <- args[1]
week            <- args[2]










### Load the corresponding imputed SNPs and SNP info for that chromosome
load(paste0("/home/phamd/imputed_snps_chr_", chr, ".RData"))
imp_snps             <- get(paste0('imp_snps_chr', chr))                    # Get imputed snp matrix
index_snpinfo        <- get(paste0('indexed_snpinfo_chr', chr))             # Get the indexed snp info dataframe
index_snpinfo$pos_bp <- round(index_snpinfo$pos * 1e6)                      # Make new column of position in bp
save_directory       <- paste0(home_directory,'week_',week,'/')             # Directory to store output









### Change directory to where fastq are store and get sample ids
setwd(fastq_directory)
sample_num     <- sapply(dir(), function(s) unlist(strsplit(s, '_'))[3])
sample_id      <- paste0('DPDP.DO2.',sample_num,'.F')







### Creating object to store results
sample_results   <- vector("list", length(dir()))
names(sample_results)   <- sample_id

                         







### Loop through each MGS_DO_* directory
for(index in 1:length(sample_num)){
  
    # Get directory name
    sample <- paste0('MGS_DO_', sample_num[index])
    sample_week_directory <- paste0(fastq_directory, sample, '/Week_', week)
  
  
  
  
    # Change to MGS_DO_* directory and week directory if it exists
    if(dir.exists(sample_week_directory)){
       setwd(sample_week_directory)
       print(sample_week_directory)
  
  
  
  
      
      
       # Read in read counts of one chromosome of one sample. Filter rows where Major and Minor allele count is 0
       sample_read_counts <- readRDS(paste0(sample,'_Week_',week,'_readcounts_chr_',chr,'.rds'))
       sample_read_counts <- sample_read_counts[!(sample_read_counts$Major_Count == 0 & sample_read_counts$Minor_Count == 0),]

  
  
  
    
      
       # Making sure all positions and snp id are the same
       aligning_pos           <- intersect(sample_read_counts$pos,index_snpinfo$pos_bp)
       filtered_index_snpinfo <- index_snpinfo[index_snpinfo$pos_bp %in% aligning_pos,]
       sample_read_counts     <- sample_read_counts[sample_read_counts$pos %in% aligning_pos,]
       imp_snps_col           <- colnames(imp_snps)[colnames(imp_snps) %in% filtered_index_snpinfo$snp]
  
  
       print(length(imp_snps_col))    
    
  
      
      
      
      
      
       # Create object to contain the results for sample pairs
       pair_results[[index]] <- array(0, dim=c(nrow(imp_snps), 3, 3, 2))
       dimnames(pair_results[[index]]) <- list(rownames(imp_snps), c("AA", "AB", "BB"),
                                            c("AA", "AB", "BB"), c("A", "B"))
  
      
      
       if(sample_id[index] %in% rownames(imp_snps)){
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
      
      } # if sample_id
    
    } # if  
    
} # loop over MB samples



saveRDS(pair_results, paste0(save_directory, "paired_results_week_",week,"_chr_", chr, ".rds"))
