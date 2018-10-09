################################################################################################################
#   This scripts takes the genotype files from /projects/churchill-lab/data/Weinstock/Pomp_Benson/genotypes to
#      remove duplicates and rename samples and marker ids. Get the phenotypes from Pomp_Benson/phenotypes to
#      check if we have genoprobs for all DO2 samples and if the sex matches.
#
#   Input:
#      1.) pomp_genoprobs_qtl2.rds
#      2.) pomp_map_qtl2.rds
#      3.) DO1DO2micecombined_Final.csv
#
#
#   Output:
#      1.) New genoprobs with sample ids and marker ids renamed
#      2.) New map with marker ids renamed
#      3.) New map as a dataframe
#
#
#   Author: Duy Pham
#   E-mail: duy.pham@jax.org
#   Date: 10/09/2018
################################################################################################################

### Load libraries
library(qtl2convert)
library(tidyverse)
library(data.table)
library(argyle)




### Options
options(stringsAsFactors = FALSE)




### Load in data
pomp_genoprobs <- readRDS('~/Desktop/Weinstock/pomp_genoprobs_qtl2.rds')
pomp_maps <- readRDS("~/Desktop/Weinstock/pomp_map_qtl2.rds")
pomp_pheno <- read.csv("~/Desktop/Weinstock/DO1DO2micecombined_Final.csv")






### Creating a new markers dataframe from the pomp_map_qtl2.rds file
#     Dimensions before: 68268 x 3
#     Dimensions after:  68268 x 5
pomp_markers <- map_list_to_df(pomp_maps)                                                    # Convert map list to markers dataframe
print(dim(pomp_markers))                                                                     # Print before dimensions
pomp_markers <- pomp_markers %>% 
                      mutate(bp = pos * 1000000, marker.id = paste0(chr,'_',bp)) %>%         # Add a bp column and marker.id column
                      select(marker.id, chr, pos, bp, marker) %>%                            # Rearranging the marker dataframe
                      dplyr::rename(orig.name = marker, marker = marker.id)                  # Reaname marker and marker.id columns
rownames(pomp_markers) <- pomp_markers$marker.id                                             # Make rownames the new marker.id column
print(dim(pomp_markers))                                                                     # Print after dimensions


### Create a new pomp_map list based on the new markers dataframe
pomp_maps <- map_df_to_list(pomp_markers, pos_column = 'pos', marker_column = 'marker')








### Checking for duplicates names. Names are generally structured as DPDP.DO[1/2].ID.sex. There are some with extended info after sex. These are
#     dupliated ones.
#
#         Length of genoprobs list: 20
#         Dimensions of genoprobs of each chromosome:  605 x 8 x #
#
#   Checking to see if all samples ID match across the genoprobs list
for(i in 2:length(pomp_genoprobs)){
  stopifnot((dimnames(pomp_genoprobs[[1]])[[1]] == dimnames(pomp_genoprobs[[i]])[[1]]))
}








### Get the names. Split the names by '.' and save as dataframe. Will use this to filter names to DO#.Sample#
orig_genoprobs_names <- dimnames(pomp_genoprobs[[1]])[[1]]
genoprobs_name <- gsub('_','.', orig_genoprobs_names)
genoprobs_name <- strsplit(genoprobs_name, '.', fixed = TRUE)
max_name_separation <- max(unlist(lapply(genoprobs_name, length)))


genoprobs_name_df <- as.data.frame(matrix(NA, nrow = length(orig_genoprobs_names), ncol = max_name_separation))
colnames(genoprobs_name_df) <- c('DPDP','DO','Sample.Number','Sex', paste0('Extra', c(1,2)))
for(i in 1:length(genoprobs_name)){
    lst <- unlist(genoprobs_name[[i]])
    genoprobs_name_df[i,1:length(lst)] <- lst
    
}







### Removing all characters after the Sex character to see duplicates. 3 Samples with duplicates. 1 in DO1 and 2 in DO2.
###   I am keeping the sample's second run: "DPDP.DO2.595.F.1" and "DPDP.DO2.420.F.1". Not sure which to keep for DO1 duplicate. Leaving it there
#     since we only have microbiomes for DO2.
filter_genoprobs_names <- gsub("*.[MF].*", "", toupper(dimnames(pomp_genoprobs[[1]])[[1]]))
dups <- c(which(duplicated(filter_genoprobs_names) | duplicated(filter_genoprobs_names, fromLast = TRUE)))   

keep = rep(TRUE,nrow(genoprobs_name_df))
for(i in seq(1,length(dups),2)){
    keep[(as.numeric(names(which.min(rowSums(is.na(genoprobs_name_df[c(dups[i],dups[i+1]),]))))))] <- FALSE
  
}







### Prepping the genoprob_name dataframe.
genoprobs_name_df <- genoprobs_name_df %>% filter(keep) %>% select(-Extra1, -Extra2) %>%
                                           mutate(Original.Name = orig_genoprobs_names[keep],
                                                  Sample.Number = as.numeric(Sample.Number),
                                                  Final.Name = paste0(DO, '.', Sample.Number),
                                                  Sex = toupper(Sex))
                              







### Change the genoprobs name: DO name and marker ID
for(i in 1:length(pomp_genoprobs)){
    pomp_genoprobs[[i]] <- pomp_genoprobs[[i]][keep,,]
    dimnames(pomp_genoprobs[[i]])[[1]] <- genoprobs_name_df$Final.Name[match(dimnames(pomp_genoprobs[[i]])[[1]], genoprobs_name_df$Original.Name)]
    dimnames(pomp_genoprobs[[i]])[[3]] <- pomp_markers$marker[match(dimnames(pomp_genoprobs[[i]])[[3]], pomp_markers$orig.name)]
}









### DO2s should be all female. All DO2s present in genoprobs. Not checking for DO1s since we don't need.
pomp_pheno <- pomp_pheno %>% filter(DO == 2)
stopifnot(sum(pomp_pheno$MouseID %in% genoprobs_name_df$Final.Name) == nrow(pomp_pheno))
stopifnot(sum(pomp_pheno$Sex == 0) == nrow(pomp_pheno))









### Save the new data
saveRDS(pomp_genoprobs, 'weinstock_DO602_genoprobs_20181009.rds')
saveRDS(pomp_markers, 'weinstock_marker_20181009.rds')
saveRDS(pomp_maps, 'weinstock_maps_20181009.rds')












### Saving the code below for future refence
# ### Predicting sex of samples that did not have same sex in genoprobs and pomp_pheno
# load("snps.megamuga.Rdata")
# geno1 <- read.beadstudio(prefix = "", snps = snps, in.path = "~/Desktop/Weinstock/UNC_PompMouse12dec2012/")
# geno2 <- read.beadstudio(prefix = "", snps = snps, in.path = "~/Desktop/Weinstock/UNC_UNL_MegaMuga31Aug2012/")
# predicted_sex <- predict.sex(cbind(geno1,geno2))
# predicted_sex <- predicted_sex[grep('DO',predicted_sex$iid),]
# predicted_sex$iid <- gsub('DP-','', predicted_sex$iid)
# predicted_sex$iid <- gsub('-','.', predicted_sex$iid)
# predicted_sex$iid <- gsub('_','.', predicted_sex$iid)
# predicted_sex$iid <- gsub('.[MF].*','', toupper(predicted_sex$iid))
# predicted_sex$iid <- sapply(predicted_sex$iid, FUN = function(x) x = paste0(strsplit(x, '.', fixed = TRUE)[[1]][1],'.',as.numeric(strsplit(x, '.', fixed = TRUE)[[1]][2])))
# predicted_sex <- predicted_sex[predicted_sex$iid %in% new_genoprobs_name,]
# predicted_sex <- predicted_sex[-which(duplicated(predicted_sex$iid)),]
# 
# pomp_pheno[pomp_pheno$MouseID %in% wrong_sex$MouseID, 'Sex'] <- 0
# sex[sex$MouseID %in% wrong_sex$MouseID,'Sex'] <- 'F'
# stopifnot(sex$Sex == genoprobs_sample_sex)

