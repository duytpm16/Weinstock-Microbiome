### Options and Libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(qtl2)
library(qtl2convert)
setwd('~/Desktop')








### Load data
genoprobs <- readRDS("~/Desktop/Weinstock Pomp Microbiome/Original Genotypes/pomp_genoprobs_qtl2.rds")
map <- readRDS("~/Desktop/Weinstock Pomp Microbiome/Original Genotypes/pomp_map_qtl2.rds")
metab <- read.delim("~/Desktop/Weinstock Pomp Microbiome/16S_taxonomic_abundance_tables_for_Pomp_DO_mouse_data/metabat_genome_cleaned_relab.tsv")
metab <- metab %>%
               remove_rownames() %>%
               column_to_rownames('genome')
markers <- map_list_to_df(map_list = map)
K <- calc_kinship(probs = genoprobs, type = 'loco', cores= 10)










### Split by Weeks
metab_w6 <- metab %>%
                  select(grep('w6', colnames(.))) %>%
                  `colnames<-`(paste0(gsub('mwgs.ain.w6.R', 'DPDP.DO2.',colnames(.), fixed = TRUE),'.F')) %>%
                  t(.)
metab_w17 <- metab %>%
                   select(grep('w17', colnames(.))) %>%
                   `colnames<-`(gsub('cholic.', '',colnames(.), fixed = TRUE)) %>%
                   `colnames<-`(gsub('protein.', '',colnames(.), fixed = TRUE)) %>%
                   `colnames<-`(paste0(gsub('mwgs.w17.R', 'DPDP.DO2.',colnames(.), fixed = TRUE),'.F')) %>%
                   t(.)
metab_w24 <- metab %>%
                   select(grep('w24', colnames(.))) %>%
                   `colnames<-`(gsub('cholic.', '',colnames(.), fixed = TRUE)) %>%
                   `colnames<-`(gsub('protein.', '',colnames(.), fixed = TRUE)) %>%
                   `colnames<-`(paste0(gsub('mwgs.w24.R', 'DPDP.DO2.',colnames(.), fixed = TRUE),'.F')) %>%
                   t(.)










# Rank Z function
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()

metab_w6_rz  <- apply(metab_w6, 2, rankZ)
metab_w17_rz <- apply(metab_w17, 2, rankZ)
metab_w24_rz <- apply(metab_w24, 2, rankZ)















### Get samples dataframe
samples <- data.frame(original_name = colnames(metab))
samples <- samples %>%
                   separate('original_name', c('mwgs','diet','week','DO_number')) %>%
                   group_by(diet, week, DO_number) %>%
                   mutate(original_name = paste(mwgs,diet,week,DO_number),
                          original_name = gsub(' ','.', original_name,fixed = TRUE),
                          mouse_id = paste0('DPDP.DO2.', gsub('R', '', DO_number, fixed = TRUE), '.F')) %>%
                   as.data.frame()

samples_w17 <- samples %>%
                       filter(grepl('w17', week)) %>%
                       select(-week) %>%
                       `rownames<-`(.$mouse_id)

samples_w24 <- samples %>%
                       filter(grepl('w24', week)) %>%
                       select(-week) %>%
                       `rownames<-`(.$mouse_id)










### Covariates for each week
covar_w17 <- model.matrix(~ diet, data = samples_w17)[,-1,drop = FALSE]
covar_w24 <- model.matrix(~ diet, data = samples_w24)[,-1,drop = FALSE]









### QTL viewer format
dataset.metab.w6 <- list(annots = data.frame(),
                         covar  = NULL,
                         covar.factors = data.frame(),
                         datatype = 'phenotype',
                         lod.peaks = list(),
                         raw = metab_w6,
                         norm = metab_w6,
                         rankz = metab_w6_rz, 
                         samples = samples_w17)
dataset.metab.w17 <- list(annots = data.frame(),
                          covar  = covar_w17,
                          covar.factors = data.frame(),
                          datatype = 'phenotype',
                          lod.peaks = list(),
                          raw = metab_w17,
                          norm = metab_w17,
                          rankz = metab_w17_rz, 
                          samples = samples_w17)
dataset.metab.w24 <- list(annots = data.frame(),
                          covar  = covar_w24,
                          covar.factors = data.frame(),
                          datatype = 'phenotype',
                          lod.peaks = list(),
                          raw = metab_w24,
                          norm = metab_w24,
                          rankz = metab_w24_rz, 
                          samples = samples_w24)








save(dataset.metab.w6, dataset.metab.w17, dataset.metab.w24, genoprobs, map , markers, K, file = 'weinstock_metab_qtl_viewer.RData')
