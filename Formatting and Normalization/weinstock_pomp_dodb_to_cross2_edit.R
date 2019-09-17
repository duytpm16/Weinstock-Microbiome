#############################################################################################
#  This script is used to edit the zip file returned by dodb.jax.org.
#        Project title: 194_Pomp_DO
#     File was downloaded on 09/14/2019. 
#     Needed to make changes to sample names so qtl2 can read in files as a cross object
#     I made a duplicated of the zip file so that I can store the original and edit the duplicate
#
#
#
#  Input:
#     map - Old map list on /projects/churchill-lab/data/Weinstock/Pomp_Benson/.
#             Shown to me by Dan Gatti.
#     covar - *_covar.csv file from dodb
#     pheno - *_pheno.csv file from dodb
#     geno  - *_geno.csv file from dodb
#
#
#  Output:
#     .Rdata - contains genoprobs in 36 state, map list, and markers dataframe 
#
#
#
#  Author: Duy Pham
#  E-mail: duy.pham@jax.org
#############################################################################################

### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(tidyverse)
library(qtl2)





### Load markers from churchill-lab HPC: /projects/churchill-lab/data/Weinstock/Pomp_Benson/genotypes/
map    <- readRDS('~/Desktop/Weinstock_Pomp_Microbiome/Genotypes/Original Genotypes/pomp_map_qtl2.rds')
marker <- map_list_to_df(map_list = map)







### Edit pheno:
#     1.) Keep samples that are DO
#     2.) Add the correct sex to sex column of covar. Previously was NAs
#     3.) Editing sample name to match sample name in geno column below
covar <- read.csv('Downloads/Pomp-194_Pomp_DO-MegaMUGA/Pomp-194_Pomp_DO-MegaMUGA_covar.csv')
covar <- covar %>% 
           filter(grepl('DP-DO[12]-', id)) %>%
           mutate(id  = gsub('Pomp|194_Pomp_DO|', '', id, fixed = TRUE),
                  sex = toupper(sapply(id, function(x) strsplit(x = x, split = '-')[[1]][4])),
                  id  = gsub('-','.',id, fixed = TRUE))
covar$sex[c(385, 579)] <- 'F'
write.csv(x = covar, file = 'Downloads/Pomp-194_Pomp_DO-MegaMUGA-2/Pomp-194_Pomp_DO-MegaMUGA_covar.csv', row.names = FALSE)









### Edit pheno:
#     1.) Keep samples that are DO
#     2.) Editing sample name to match sample name in geno column below
pheno <- read.csv('~/Downloads/Pomp-194_Pomp_DO-MegaMUGA/Pomp-194_Pomp_DO-MegaMUGA_pheno.csv')
pheno <- pheno %>%
           filter(grepl('DP-DO[12]-', id)) %>%
           mutate(id  = gsub('Pomp|194_Pomp_DO|', '', id, fixed = TRUE),
                  id  = gsub('-','.',id, fixed = TRUE))  
write.csv(x = covar, file = '~/Downloads/Pomp-194_Pomp_DO-MegaMUGA-2/Pomp-194_Pomp_DO-MegaMUGA_pheno.csv', row.names = FALSE)








### Edit genoprobs:
#     1.) Keep columns that are marker and DO samples
#     2.) Keep markers that are in original markers dataframe
geno <- read.csv('~/Downloads/Pomp-194_Pomp_DO-MegaMUGA/Pomp-194_Pomp_DO-MegaMUGA_geno.csv')
geno <- geno[,c('marker', grep('.DO[12].', colnames(geno), value = TRUE))]
geno <- geno[geno$marker %in% marker$marker,]
write.csv(x = geno, 'Downloads/Pomp-194_Pomp_DO-MegaMUGA-2/Pomp-194_Pomp_DO-MegaMUGA_geno.csv', row.names = FALSE)







### Create cross and genoprobs object
cross    <- read_cross2(file = '~/Downloads/Pomp-194_Pomp_DO-MegaMUGA-2/Pomp-194_Pomp_DO-MegaMUGA.json')
genoprob <- calc_genoprob(cross = cross, map = map, cores = 0)





### Save
save(map, marker, genoprob, file = '~/Desktop/weinstock_pomp_36_state_genoprob.Rdata')
