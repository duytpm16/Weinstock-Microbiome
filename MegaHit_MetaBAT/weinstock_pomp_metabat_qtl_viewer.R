### Options and Libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(qtl2)
library(qtl2convert)
setwd('~/Desktop')








### Load data
metadata <- read.csv('~/Desktop/Weinstock Pomp Microbiome/Phenotypes/DO1DO2micecombined_Final.csv')
taxonomy <- read.delim('~/Desktop/Weinstock Pomp Microbiome/MetaBAT/megahit_metabat_taxonomy.tsv', sep = '|', header = FALSE)
metabat <- read.delim("~/Desktop/Weinstock Pomp Microbiome/MetaBAT/metabat_genome_cleaned_relab.tsv")
metabat <- metabat %>%
               remove_rownames() %>%
               column_to_rownames('genome') %>%
               `rownames<-`(gsub('.','_',rownames(.), fixed = TRUE))











### Split by Weeks
metabat_w6 <- metabat %>%
                  select(grep('w6', colnames(.))) %>%
                  `colnames<-`(paste0(gsub('mwgs.ain.w6.R', 'DO2.',colnames(.), fixed = TRUE))) %>%
                  `colnames<-`(gsub('.', '-', colnames(.), fixed = TRUE)) %>%
                  t(.)

metabat_w17 <- metabat %>%
                   select(grep('w17', colnames(.))) %>%
                   `colnames<-`(gsub('cholic.', '',colnames(.), fixed = TRUE)) %>%
                   `colnames<-`(gsub('protein.', '',colnames(.), fixed = TRUE)) %>%
                   `colnames<-`(paste0(gsub('mwgs.w17.R', 'DO2.',colnames(.), fixed = TRUE))) %>%
                   `colnames<-`(gsub('.', '-', colnames(.), fixed = TRUE)) %>%
                   t(.)

metabat_w24 <- metabat %>%
                   select(grep('w24', colnames(.))) %>%
                   `colnames<-`(gsub('cholic.', '',colnames(.), fixed = TRUE)) %>%
                   `colnames<-`(gsub('protein.', '',colnames(.), fixed = TRUE)) %>%
                   `colnames<-`(paste0(gsub('mwgs.w24.R', 'DO2.',colnames(.), fixed = TRUE))) %>%
                   `colnames<-`(gsub('.', '-', colnames(.), fixed = TRUE)) %>%
                   t(.)

intersect_samples   <- intersect(rownames(metabat_w6), rownames(metabat_w17))
intersect_samples   <- intersect(intersect_samples, rownames(metabat_w24))

metabat_w6  <- metabat_w6[intersect_samples,]
metabat_w17 <- metabat_w17[intersect_samples,]
metabat_w24 <- metabat_w24[intersect_samples,]





### Change
metabat_w6_w17  <- (metabat_w6 - metabat_w17) 
metabat_w17_w24 <- (metabat_w17 - metabat_w24) 






# Rank Z function
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

metabat_w6_rz   <- apply(metabat_w6,  2, rankZ)
metabat_w17_rz  <- apply(metabat_w17, 2, rankZ)
metabat_w24_rz  <- apply(metabat_w24, 2, rankZ)
metabat_w6_w17_rz  <- apply(metabat_w6_w17, 2, rankZ)
metabat_w17_w24_rz <- apply(metabat_w17_w24,2, rankZ)














### Fix Taxonomy Table
colnames(taxonomy) <- c('MetaBat.ID','Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
taxonomy$MetaBat.ID <- gsub('.', '_', sapply(taxonomy$MetaBat.ID, FUN = function(x) strsplit(x, split = '\t', fixed = TRUE)[[1]][1]), fixed = TRUE)

taxonomy[16,] <- c(taxonomy[16, -2], taxonomy[17,1])
taxonomy      <- taxonomy[-17,]
taxonomy      <- taxonomy[taxonomy$MetaBat.ID %in% rownames(metabat),]

bin_order <- data.frame(index = 1:nrow(taxonomy), id = taxonomy$MetaBat.ID) 
bin_order <- bin_order %>%
                       tidyr::separate(id, into = c('megahit','metebat','bin','number'), sep = '_') %>%
                       mutate(number = as.numeric(number)) %>%
                       dplyr:: arrange(number)


taxonomy <- taxonomy[bin_order$index, ]




for(i in 1:nrow(taxonomy)){

    p <- if(length(grep('p__', taxonomy[i,2:7], value = TRUE)) != 0) grep('p__', taxonomy[i,2:7], value = TRUE) else NA
    c <- if(length(grep('c__', taxonomy[i,2:7], value = TRUE)) != 0) grep('c__', taxonomy[i,2:7], value = TRUE) else NA
    o <- if(length(grep('o__', taxonomy[i,2:7], value = TRUE)) != 0) grep('o__', taxonomy[i,2:7], value = TRUE) else NA
    f <- if(length(grep('f__', taxonomy[i,2:7], value = TRUE)) != 0) grep('f__', taxonomy[i,2:7], value = TRUE) else NA
    g <- if(length(grep('g__', taxonomy[i,2:7], value = TRUE)) != 0) grep('g__', taxonomy[i,2:7], value = TRUE) else NA
    s <- if(length(grep('s__', taxonomy[i,2:7], value = TRUE)) != 0) grep('s__', taxonomy[i,2:7], value = TRUE) else NA


    taxonomy[i, 'Phylum']  <- gsub("[[:punct:]]", '', gsub("_", ' ', gsub('p__', '', p, fixed = TRUE), fixed = TRUE))
    taxonomy[i, 'Class']   <- gsub("[[:punct:]]", '', gsub("_", ' ', gsub('c__', '', c, fixed = TRUE), fixed = TRUE))
    taxonomy[i, 'Order']   <- gsub("[[:punct:]]", '', gsub("_", ' ', gsub('o__', '', o, fixed = TRUE), fixed = TRUE))
    taxonomy[i, 'Family']  <- gsub("[[:punct:]]", '', gsub("_", ' ', gsub('f__', '', f, fixed = TRUE), fixed = TRUE))
    taxonomy[i, 'Genus']   <- gsub("[[:punct:]]", '', gsub("_", ' ', gsub('g__', '', g, fixed = TRUE), fixed = TRUE))
    taxonomy[i, 'Species'] <- gsub("[[:punct:]]", '', gsub("_", ' ', gsub('s__', '', s, fixed = TRUE), fixed = TRUE))
}













### Make samples dataframe
metadata$MouseID <- gsub('.','-',metadata$MouseID, fixed = TRUE)
metadata$CtClr[metadata$CtClr == 1] <- 'White'
metadata$CtClr[metadata$CtClr == 2] <- 'Agouti'
metadata$CtClr[metadata$CtClr == 3] <- 'Black'
metadata$CtClr[metadata$CtClr == 4] <- 'Other'

samples <- data.frame(original_name = colnames(metabat))
samples <- samples %>%
                   separate('original_name', c('mwgs','diet','week','DO_number')) %>%
                   group_by(diet, week, DO_number) %>%
                   mutate(original_name = paste(mwgs,diet,week,DO_number),
                          original.name = gsub(' ','.', original_name),
                          mouse.id = paste0('DO2-', gsub('R', '', DO_number, fixed = TRUE)),
                          sex = 'F',
                          generation = 11,
                          coat.color = metadata$CtClr[match(mouse.id, metadata$MouseID)],
                          age.death = metadata$AgeSac[match(mouse.id, metadata$MouseID)]) %>%
                   as.data.frame()
samples <- samples[samples$mouse.id %in% intersect_samples,]

samples_w6 <- samples %>%
                      filter(grepl('w6', week)) %>%
                      select(mouse.id, sex, diet, coat.color, generation, age.death, original.name) %>%
                      `rownames<-`(.$mouse.id)

samples_w17 <- samples %>%
                       filter(grepl('w17', week)) %>%
                       select(mouse.id, sex, diet, coat.color, generation, age.death, original.name) %>%
                       `rownames<-`(.$mouse.id)

samples_w24 <- samples %>%
                       filter(grepl('w24', week)) %>%
                       select(mouse.id, sex, diet, coat.color, generation, age.death, original.name) %>%
                       `rownames<-`(.$mouse.id)













### Covariates for each week
covar_w17 <- model.matrix(~ diet, data = samples_w17)[intersect_samples,-1,drop = FALSE]
covar_w24 <- model.matrix(~ diet, data = samples_w24)[intersect_samples,-1,drop = FALSE]


covar.factors <- data.frame(column.name = c('diet'),
                            diplay.name = c('Diet'),
                            int.covar   = c('factor'),
                            lod.peaks   = c('diet_int'),
                            covar.name  = c('diet'))













### QTL viewer format
dataset.metabat.w6 <- list(annots      = taxonomy,
                         covar         = NULL,
                         covar.factors = NULL,
                         datatype      = 'phenotype',
                         lod.peaks     = list(),
                         raw           = metabat_w6[,taxonomy$MetaBat.ID],
                         norm          = metabat_w6[,taxonomy$MetaBat.ID],
                         rankz         = metabat_w6_rz[,taxonomy$MetaBat.ID],
                         samples       = samples_w17)

dataset.metabat.w17 <- list(annots      = taxonomy,
                          covar         = covar_w17,
                          covar.factors = covar.factors,
                          datatype      = 'phenotype',
                          lod.peaks     = list(),
                          raw           = metabat_w17[,taxonomy$MetaBat.ID],
                          norm          = metabat_w17[,taxonomy$MetaBat.ID],
                          rankz         = metabat_w17_rz[,taxonomy$MetaBat.ID],
                          samples       = samples_w17)

dataset.metabat.w24 <- list(annots      = taxonomy,
                          covar         = covar_w24,
                          covar.factors = covar.factors,
                          datatype      = 'phenotype',
                          lod.peaks     = list(),
                          raw           = metabat_w24[,taxonomy$MetaBat.ID],
                          norm          = metabat_w24[,taxonomy$MetaBat.ID],
                          rankz         = metabat_w24_rz[,taxonomy$MetaBat.ID],
                          samples       = samples_w24)



dataset.metabat.w6.w17.change <- list(annots        = taxonomy,
                                      covar         = covar_w17,
                                      covar.factors = covar.factors,
                                      datatype      = 'phenotype',
                                      lod.peaks     = list(),
                                      raw.w6        = metabat_w6[,taxonomy$MetaBat.ID],
                                      raw.w17       = metabat_w17[,taxonomy$MetaBat.ID],
                                      change        = metabat_w6_w17[,taxonomy$MetaBat.ID],
                                      rankz         = metabat_w6_w17_rz[,taxonomy$MetaBat.ID],
                                      samples       = samples_w17)

dataset.metabat.w17.w24.change <- list(annots        = taxonomy,
                                       covar         = covar_w24,
                                       covar.factors = covar.factors,
                                       datatype      = 'phenotype',
                                       lod.peaks     = list(),
                                       raw.w17       = metabat_w17[,taxonomy$MetaBat.ID],
                                       raw.w24       = metabat_w24[,taxonomy$MetaBat.ID],
                                       change        = metabat_w17_w24[,taxonomy$MetaBat.ID],
                                       rankz         = metabat_w17_w24_rz[,taxonomy$MetaBat.ID],
                                       samples       = samples_w24)




rm(list = ls()[!grepl('dataset.', ls())])
save.image('~/Desktop/Weinstock Pomp Microbiome/MetaBat/weinstock_metabat_qtl_viewer.RData')
