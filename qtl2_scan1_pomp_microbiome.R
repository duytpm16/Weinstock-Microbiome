### This script reads in the required input data file generated from gather_qtl2_scan1_input_data.R
#       to run the qtl2 scan1 function.
#   The QTL scan can be ran in 'chunks' or all at once.
#       Example for chunk size/number: 
#           Suppose there are 5433 phenotype columns. If chunk size is 1000, then there should be 6 different chunks
#             of scan1 runs, with the chunk_number value being 1-6 to get the column numbers:
#             1-1000,1001-2000,2001-3000,3001-4000,4001-5000,5001-5433, respectively.
#
#           If you do not want to run in chunks, set use_chunks to FALSE.
#
### Input:
#       1: input.file:    Path + prefix to the qtl2 input data generated from gather_qtl2_scan1_input_data.R
#       2: num_cores:     Number of cores to run
#       3: should_rankz:  Logical value to use the rankz dataset instead of normalized
#       4: use_chunks:    Logical value to run QTL scans in chunks
#       5: use_int:       Logical value to use an interaction term
#       6: chunk_number:  Numeric value of the chunk number. Not needed if use_chunks is FALSE
#       7: chunk_size:    Numeric value of chunk size. Should be consistent. Not needed if use_chunks is FALSE
#       8: int_name:      Name of the interaction term. Not needed if use_int is FALSE
#
### Output: 
#       1: Matrix containing LOD scoress for each of the phenotype that was given to scan1 at each marker.
#
### Author: Duy Pham, 'phenotype range run' was taken from Dan Gatti
### Date:   July 10, 2018
### E-mail: duy.pham@jax.org
####################################################################################################################



### Install required library packages
# install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl", "dplyr))
# library(devtools)
# install_github("rqtl/qtl2", "rqtl/qtl2convert")

options(stringsAsFactors = FALSE)

### Load required library packages
library(qtl2) 


### Command line arguments / Variables to change
# 1: input.file:    Path + prefix to the qtl2 input data generated from gather_qtl2_scan1_input_data.R
# 2: num_cores:     Number of cores to run
# 3: should_rankz:  Logical value to use the rankz dataset instead of normalized
# 4: use_chunks:    Logical value to run QTL scans in chunks
# 5: use_int:       Logical value to use an interaction term
# 6: chunk_number:  Numeric value of the chunk number. Not needed if use_chunks is FALSE
# 8: int_name:      Name of the interaction term. Not needed if use_int is FALSE
args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.RData"))
taxa <- args[3]
num_cores <- as.numeric(args[2])
should_rankz <- as.logical(args[3])



### Check to see if required data are loaded
stopifnot(c("genoprobs", "K", "map", "markers") %in% ls())



### If should_rankz is true used the normalized rankz dataset, else use the normalized dataset.
if(should_rankz){
   w6 <- get(paste0('rZ_',taxa,'_w6'))
   w17 <- get(paste0('rZ_',taxa,'_w17'))
   w24 <- get(paste0('rZ_',taxa,'_w24'))
  
}else{
   w6 <- get(paste0(taxa,'_w6'))
   w17 <- get(paste0(taxa,'_w17'))
   w24 <- get(paste0(taxa,'_w24'))
}






### Running QTL2 scan1 function


qtl_w6 <- scan1(genoprobs = genoprobs, 
                pheno = w6,
                kinship = K, 
                addcovar = covar_w6,
                intcovar = NULL, 
                cores = num_cores)

qtl_w17 <- scan1(genoprobs = genoprobs, 
                pheno = w17,
                kinship = K, 
                addcovar = covar_w17,
                intcovar = NULL, 
                cores = num_cores)

qtl_w24 <- scan1(genoprobs = genoprobs, 
                pheno = w24,
                kinship = K, 
                addcovar = covar_w,
                intcovar = NULL, 
                cores = num_cores)



assign(x = paste0(taxa,'_qtl'),
       value = list(w6 = qtl_w6, w17 = qtl_w17, w24 = qtl_w24), 
       envir = .GlobalEnv)



  
# Save output of scan1 to current directory
if(should_rankz){
   saveRDS(get(paste0(taxa,'_qtl')), file = paste0(args[1], "_rZ_qtl_lod.rds"))
}else{
   saveRDS(get(paste0(taxa,'_qtl')), file = paste0(args[1], "_norm_qtl_lod.rds"))
}
