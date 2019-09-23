#################################################################################
#   This script is used to modify the genoprobs and map original .rds file.
#     1.) I am changing the sample name to DO2.# instead of DPDP.DO2.#.F
#     2.) Remove all DO1 cohort samples to reduce memory. There were no microbiome counts for these samples
#     3.) Create a markers dataframe (map list -> dataframe)
#     4.) Calculate kinship matrix
#
#   Genoprobs obtained on cadillac HPC: /projects/churchill-lab/data/Weinstock/Pomp_Benson/genotypes/pomp_genoprobs_qtl2.rds
#   map list obtained on cadillac HPC: /projects/churchill-lab/data/Weinstock/Pomp_Benson/genotypes/pomp_map_qtl2.rds
#
#
#
#   Input:
#     genoprobs - obtained from directory above
#     map - obtained from directory above
#
#   Outut:
#     .Rdata file containing genoprobs, map, marker dataframe, and kinship list for qtl2
#
#
#
#   Author: Duy Pham
#   E-mail: duy.pham@jax.org
#################################################################################

### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(tidyverse)
library(qtl2)





### Load in original data
genoprobs <- readRDS('~/Desktop/Weinstock_Pomp_Microbiome/Genotypes/Original Genotypes/pomp_genoprobs_qtl2.rds')
map <- readRDS('~/Desktop/Weinstock_Pomp_Microbiome/Genotypes/Original Genotypes/pomp_map_qtl2.rds')







### Editing genoprobs
#    1.) Remove DO1 cohort to reduce memory
#    2.) Edit sample names to DO2.# format
genoprobs <- probs_qtl2_to_doqtl(probs = genoprobs)
genoprobs <- genoprobs[grep('DO2', dimnames(genoprobs)[[1]]),,]
dimnames(genoprobs)[[1]] <- gsub('DPDP[.]|[.]F', '', dimnames(genoprobs)[[1]])






### Create a maps dataframe
markers <- map_list_to_df(map_list = map, chr_column = 'chr', pos_column = 'pos', marker_column = 'marker.id')
markers <- markers %>% select(marker.id, chr, pos)









### Put genoprobs back to qtl2 format and create a kinship matrix
genoprobs <- probs_doqtl_to_qtl2(probs = genoprobs, map = markers, 
                                 chr_column = 'chr', pos_column = 'pos', marker_column = 'marker.id')

K <- calc_kinship(probs = genoprobs, type = 'loco', cores = 0)











### Save
save(genoprobs, K, map, markers, file = 'weinstock_pomp_original_genoprobs.Rdata')
