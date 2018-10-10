### This script reads in the required input data file generated from gather_microbiome_qtl2_input.R
#       to run the qtl2 scan1 function.
#
#
### Input:
#       1: input.file:    Path + prefix to the qtl2 input data generated from gather_microbiome_qtl2_input.R
#       2: taxa:          which taxa to perform qtl2 scan on
#       3: num_cores:     Number of cores to run
#       4: should_rankz:  Logical value to use the rankz dataset instead of normalized

#
### Output: 
#       1: List of LOD matrices for w6, w17, and w24 for the given taxa in Input (2).
#
### Author: Duy Pham
### Date:   October 10, 2018
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
# 2: taxa:          which taxa to perform qtl2 scan on
# 3: num_cores:     Number of cores to run
# 4: should_rankz:  Logical value to use the rankz dataset instead of normalized

args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.RData"))
taxa <- args[2]
num_cores <- as.numeric(args[3])
should_rankz <- as.logical(args[4])



### Check to see if required data are loaded
stopifnot(c("genoprobs", "K", "map", "markers") %in% ls())



### If should_rankz is true used the normalized rankz dataset, else use the normalized dataset.
if(should_rankz){
   w6 <- paste0('rZ_',taxa,'_w6')
   w17 <- paste0('rZ_',taxa,'_w17')
   w24 <- paste0('rZ_',taxa,'_w24')
  
}else{
   w6 <- paste0(taxa,'_w6')
   w17 <- paste0(taxa,'_w17')
   w24 <- paste0(taxa,'_w24')
}






### Running QTL2 scan1 function


qtl_w6 <- scan1(genoprobs = genoprobs, 
                pheno = expr,
                kinship = K, 
                addcovar = covar_w6,
                intcovar = NULL, 
                cores = num_cores)

qtl_w17 <- scan1(genoprobs = genoprobs, 
                pheno = expr,
                kinship = K, 
                addcovar = covar_w17,
                intcovar = NULL, 
                cores = num_cores)

qtl_w24 <- scan1(genoprobs = genoprobs, 
                pheno = expr,
                kinship = K, 
                addcovar = covar_w6,
                intcovar = NULL, 
                cores = num_cores)



assign(x = paste0(taxa,'_qtl'),
       value = list(paste0(taxa,'_w6') = qtl_w6, paste0(taxa,'_w17') = qtl_w17, paste0(taxa,'_w24') = qtl_w24), 
       envir = .GlobalEnv)



  
# Save output of scan1 to current directory
if(should_rankz){
   saveRDS(get(paste0(taxa,'_qtl')), file = paste0(args[1],'_',taxa, "_rZ_qtl_lod.rds"))
}else{
   saveRDS(get(paste0(taxa,'_qtl')), file = paste0(args[1], '_',taxa, "_norm_qtl_lod.rds"))
}
