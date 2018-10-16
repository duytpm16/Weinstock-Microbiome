
### Load library packages
options(stringsAsFactors = FALSE)

library(openxlsx)
library(dplyr)



### Read in the microbiome data that have been normalized by the median-normalized method.
genus_16s <- read.xlsx('~/Desktop/Weinstock Microbiome/16S_taxonomic_abundance_tables_for_Pomp_DO_mouse_data/taxonomic_abundance_filtered_mediannormalized_20181004.xlsx',
                       sheet = 'genus_normalized_abundance')
family_16s <- read.xlsx('~/Desktop/Weinstock Microbiome/16S_taxonomic_abundance_tables_for_Pomp_DO_mouse_data/taxonomic_abundance_filtered_mediannormalized_20181004.xlsx',
                        sheet = 'family_normalized_abundance')
order_16s <- read.xlsx('~/Desktop/Weinstock Microbiome/16S_taxonomic_abundance_tables_for_Pomp_DO_mouse_data/taxonomic_abundance_filtered_mediannormalized_20181004.xlsx',
                       sheet = 'order_normalized_abundance')
class_16s <- read.xlsx('~/Desktop/Weinstock Microbiome/16S_taxonomic_abundance_tables_for_Pomp_DO_mouse_data/taxonomic_abundance_filtered_mediannormalized_20181004.xlsx',
                       sheet = 'class_normalized_abundance')
phylum_16s <- read.xlsx('~/Desktop/Weinstock Microbiome/16S_taxonomic_abundance_tables_for_Pomp_DO_mouse_data/taxonomic_abundance_filtered_mediannormalized_20181004.xlsx',
                        sheet = 'phylum_normalized_abundance')

samples_16s <- read.xlsx('~/Desktop/Weinstock Microbiome/16S_taxonomic_abundance_tables_for_Pomp_DO_mouse_data/taxonomic_abundance_filtered_mediannormalized_20181004.xlsx',
                         sheet = 'Sample IDs')









### Change Mouse name and reorder dataframe
samples_16s <- samples_16s %>%
                    filter(!duplicated(.)) %>%
                    mutate(Mouse = paste0('DO2.',Mouse)) %>%
                    arrange(Mouse, SampleName, Diet)








### Filter normalized data by week and rename 
samples_w6 <- samples_16s %>% filter(Week == 6) %>% `rownames<-`(.$Mouse)
samples_w17 <- samples_16s %>% filter(Week == 17) %>% `rownames<-`(.$Mouse)
samples_w24 <- samples_16s %>% filter(Week == 24) %>% `rownames<-`(.$Mouse)



genus_16s <- genus_16s %>% `rownames<-`(.[,'Taxon']) %>% select(-Taxon) %>% t() %>% `rownames<-`(gsub('16s-', '', rownames(.)))
genus_w6 <- genus_16s[grep('w6',rownames(genus_16s)),] %>% `rownames<-`(samples_w6$Mouse[match(samples_w6$SampleName, rownames(.))])
genus_w17 <- genus_16s[grep('w17',rownames(genus_16s)),] %>% `rownames<-`(samples_w17$Mouse[match(samples_w17$SampleName, rownames(.))])
genus_w24 <- genus_16s[grep('w24',rownames(genus_16s)),] %>% `rownames<-`(samples_w24$Mouse[match(samples_w24$SampleName, rownames(.))])


family_16s <- family_16s %>% `rownames<-`(.[,'Taxon']) %>% select(-Taxon) %>% t() %>% `rownames<-`(gsub('16s-', '', rownames(.)))
family_w6 <- family_16s[grep('w6',rownames(family_16s)),] %>% `rownames<-`(samples_w6$Mouse[match(samples_w6$SampleName, rownames(.))])
family_w17 <- family_16s[grep('w17',rownames(family_16s)),] %>% `rownames<-`(samples_w17$Mouse[match(samples_w17$SampleName, rownames(.))])
family_w24 <- family_16s[grep('w24',rownames(family_16s)),] %>% `rownames<-`(samples_w24$Mouse[match(samples_w24$SampleName, rownames(.))])


order_16s <- order_16s %>% `rownames<-`(.[,'Taxon']) %>% select(-Taxon) %>% t() %>% `rownames<-`(gsub('16s-', '', rownames(.)))
order_w6 <- order_16s[grep('w6',rownames(order_16s)),] %>% `rownames<-`(samples_w6$Mouse[match(samples_w6$SampleName, rownames(.))])
order_w17 <- order_16s[grep('w17',rownames(order_16s)),] %>% `rownames<-`(samples_w17$Mouse[match(samples_w17$SampleName, rownames(.))])
order_w24 <- order_16s[grep('w24',rownames(order_16s)),] %>% `rownames<-`(samples_w24$Mouse[match(samples_w24$SampleName, rownames(.))])


class_16s <- class_16s %>% `rownames<-`(.[,'Taxon']) %>% select(-Taxon) %>% t() %>% `rownames<-`(gsub('16s-', '', rownames(.)))
class_w6 <- class_16s[grep('w6',rownames(class_16s)),] %>% `rownames<-`(samples_w6$Mouse[match(samples_w6$SampleName, rownames(.))])
class_w17 <- class_16s[grep('w17',rownames(class_16s)),] %>% `rownames<-`(samples_w17$Mouse[match(samples_w17$SampleName, rownames(.))])
class_w24 <- class_16s[grep('w24',rownames(class_16s)),] %>% `rownames<-`(samples_w24$Mouse[match(samples_w24$SampleName, rownames(.))])


phylum_16s <- phylum_16s %>% `rownames<-`(.[,'Taxon']) %>% select(-Taxon) %>% t() %>% `rownames<-`(gsub('16s-', '', rownames(.)))
phylum_w6 <- phylum_16s[grep('w6',rownames(phylum_16s)),] %>% `rownames<-`(samples_w6$Mouse[match(samples_w6$SampleName, rownames(.))])
phylum_w17 <- phylum_16s[grep('w17',rownames(phylum_16s)),] %>% `rownames<-`(samples_w17$Mouse[match(samples_w17$SampleName, rownames(.))])
phylum_w24 <- phylum_16s[grep('w24',rownames(phylum_16s)),] %>% `rownames<-`(samples_w24$Mouse[match(samples_w24$SampleName, rownames(.))])




### Create covar
covar_w6 <- NULL
covar_w17 <- model.matrix(~Diet, samples_w17)[,-1, drop = FALSE]
covar_w24 <- model.matrix(~Diet, samples_w24)[,-1, drop = FALSE]






## Making sure all names are in both covar and taxa dataset
for(i in c(17,24)){
    stopifnot(sum(rownames(get(paste0('covar_w',i))) %in% rownames(get(paste0('genus_w',i)))) == nrow(get(paste0('genus_w',i))))
    stopifnot(sum(rownames(get(paste0('covar_w',i))) %in% rownames(get(paste0('family_w',i)))) == nrow(get(paste0('family_w',i))))
    stopifnot(sum(rownames(get(paste0('covar_w',i))) %in% rownames(get(paste0('order_w',i)))) == nrow(get(paste0('order_w',i))))
    stopifnot(sum(rownames(get(paste0('covar_w',i))) %in% rownames(get(paste0('class_w',i)))) == nrow(get(paste0('class_w',i))))
    stopifnot(sum(rownames(get(paste0('covar_w',i))) %in% rownames(get(paste0('phylum_w',i)))) == nrow(get(paste0('phylum_w',i))))
}







### Read in genoprobs, markers, and map data. Calculate Kinship.
genoprobs <- readRDS("~/Desktop/Weinstock Microbiome/Filter Genotypes/weinstock_DO602_genoprobs_20181009.rds")
markers <- readRDS("~/Desktop/Weinstock Microbiome/Filter Genotypes/weinstock_markers_20181009.rds")
map <- readRDS("~/Desktop/Weinstock Microbiome/Filter Genotypes/weinstock_maps_20181009.rds")
K <- calc_kinship(genoprobs, type = 'loco', cores = 0)







### RankZ the data
rankZ = function(x) {
      x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
      return(qnorm(x))
} # rankZ()


rZ_genus_w6 <- apply(genus_w6, 2, rankZ)
rZ_genus_w17 <- apply(genus_w17, 2, rankZ)
rZ_genus_w24 <- apply(genus_w24, 2, rankZ)


rZ_family_w6 <- apply(family_w6, 2, rankZ)
rZ_family_w17 <- apply(family_w17, 2, rankZ)
rZ_family_w24 <- apply(family_w24, 2, rankZ)


rZ_order_w6 <- apply(order_w6, 2, rankZ)
rZ_order_w17 <- apply(order_w17, 2, rankZ)
rZ_order_w24 <-  apply(order_w24, 2, rankZ)


rZ_class_w6 <- apply(class_w6, 2, rankZ)
rZ_class_w17 <- apply(class_w17, 2, rankZ)
rZ_class_w24 <- apply(class_w24, 2, rankZ)


rZ_phylum_w6 <- apply(phylum_w6, 2, rankZ)
rZ_phylum_w17 <- apply(phylum_w17, 2, rankZ)
rZ_phylum_w24 <- apply(phylum_w24, 2, rankZ)


rm(class_16s, family_16s, genus_16s, order_16s, phylum_16s, samples_16s, rankZ,i)


save.image(file = 'weinstock_16s_microbiome_qtl2_input.RData')
