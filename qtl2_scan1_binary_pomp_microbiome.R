### This script reads in the required input data file generated from gather_microbiome_qtl2_input.R
#       which contains the microbiome abundancies as binary matrices
#
#
### Input:
#       1: input.file:    Path + prefix to the qtl2 input data generated from gather_microbiome_qtl2_input.R
#       2: taxa:          Which taxa to run qtl scan on
#       3: num_cores:     Number of cores to run
#       4: should_rankz:  Logical value to use the rankz dataset instead of normalized

#
### Output: 
#       1: List containing matrices of LOD scores for each week of a given taxa.
#
### Author: Duy Pham
### Date:   October 10, 2018
### E-mail: duy.pham@jax.org
####################################################################################################################


### Install required library packages
# install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl", "dplyr))
# library(devtools)
# install_github("rqtl/qtl2")

options(stringsAsFactors = FALSE)

### Load required library packages
library(qtl2) 







### Command line arguments / Variables to change
# 1: input.file:    Path + prefix to the qtl2 input data generated from gather_microbiome_qtl2_input.R
# 2: taxa:          Which taxa to run qtl2 scan1 function on 
# 3: num_cores:     Number of cores to run
# 4: should_rankz:  Logical value to use the rankz dataset instead of normalized

args = commandArgs(trailingOnly = TRUE)

load(paste0(args[1],"_qtl2_input.RData"))
taxa <- args[2]
num_cores <- as.numeric(args[3])
should_rankz <- as.logical(args[4])








### Check to see if required data are loaded
stopifnot(c("genoprobs", "K", "map", "markers") %in% ls())







### Read in binary matrix of given taxa
w6 <- get(paste0(taxa,'_w6'))
w17 <- get(paste0(taxa,'_w17'))
w24 <- get(paste0(taxa,'_w24'))








### Running QTL2 scan1 function
qtl_w6 <- scan1(genoprobs = genoprobs, 
                pheno = w6,
                kinship = K, 
                addcovar = covar_w6,
                intcovar = NULL,
                model = "binary",
                cores = num_cores)

qtl_w17 <- scan1(genoprobs = genoprobs, 
                 pheno = w17,
                 kinship = K, 
                 addcovar = covar_w17,
                 intcovar = NULL,
                 model = "binary",
                 cores = num_cores)

qtl_w24 <- scan1(genoprobs = genoprobs, 
                 pheno = w24,
                 kinship = K, 
                 addcovar = covar_w24,
                 intcovar = NULL,
                 model = "binary",
                 cores = num_cores)


qtl_list <- list(w6 = qtl_w6, w17 = qtl_w17, w24 = qtl_w24)







# Save output of scan1 to current directory
saveRDS(qtl_list, file = paste0(args[1], '_',taxa,"_binary_qtl_lod.rds"))
