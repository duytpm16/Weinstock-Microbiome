### Options and Libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(qtl2convert)
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(qtl2)
select <- dplyr::select







### Load data
load("~/Desktop/Weinstock_Pomp_Microbiome/Genotypes/Original Genotypes/weinstock_pomp_genoprobs_map_markers_K.RData")
metadata <- read.csv("~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/Mouse Info/DO1DO2micecombined_Final.csv")
sampleSheet <- read.xlsx(xlsxFile = "~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/16S/Raw Counts/taxonomic_abundance_filtered_unnormalized_20181004.xlsx")
taxa_table  <- read.xlsx(xlsxFile = "~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/16S/Raw Counts/taxonomic_abundance_filtered_unnormalized_20181004.xlsx", sheet = 'OTU Taxonomy')
otu_count   <- read.xlsx(xlsxFile = "~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/16S/Raw Counts/taxonomic_abundance_filtered_unnormalized_20181004.xlsx", sheet = 'otu_abundance')









### Editing taxa table
colnames(taxa_table) <- c('OTU', 'domain', 'phylum', 'class', 'order', 'family', 'genus')
taxa_table <- apply(taxa_table, 2, FUN = function(x) gsub('"', "",  x, fixed = TRUE))
taxa_table <- apply(taxa_table, 2, FUN = function(x) gsub('_', " ", x, fixed = TRUE))
taxa_table <- taxa_table %>% as.data.frame() %>% select(OTU, genus, family, order, class, phylum, domain)
taxa_table$OTU <- paste0('OTU-', gsub('usearchOTU','',taxa_table$OTU))







### Adding coat color to metadata
metadata$MouseID <- gsub('.','-',metadata$MouseID, fixed = TRUE)
metadata$CtClr[metadata$CtClr == 1] <- 'white'
metadata$CtClr[metadata$CtClr == 2] <- 'agouti'
metadata$CtClr[metadata$CtClr == 3] <- 'black'
metadata$CtClr[metadata$CtClr == 4] <- 'other'








### Create an annot.sample for each week
sampleSheet <- sampleSheet[!duplicated(sampleSheet),]  %>% 
                  arrange(Mouse, Week) %>% 
                  mutate(mouse.id   = paste0('DO2-', Mouse), 
                         sex        = 'F', 
                         generation = '11',
                         coat.color = metadata$CtClr[match(mouse.id, metadata$MouseID)],
                         age.of.death = metadata$AgeSac[match(mouse.id, metadata$MouseID)],
                         original.id = paste0('16s-', SampleName)) %>% 
                  dplyr::rename(diet = Diet,
                                week = Week) %>%
                  select(mouse.id, coat.color, sex, diet, week, generation, age.of.death, original.id)

sampleSheet_w6  <- sampleSheet[sampleSheet$week == 6,] %>% mutate(diet = factor(diet))
sampleSheet_w17 <- sampleSheet[sampleSheet$week == 17,] %>% mutate(diet = factor(diet))
sampleSheet_w24 <- sampleSheet[sampleSheet$week == 24,] %>% mutate(diet = factor(diet))

intersect_samples <- intersect(intersect(sampleSheet_w6$mouse.id, sampleSheet_w17$mouse.id), sampleSheet_w24$mouse.id)
intersect_samples <- intersect_samples[intersect_samples %in% dimnames(genoprobs[[1]])[[1]]]

sampleSheet_w6  <- sampleSheet_w6[sampleSheet_w6$mouse.id %in% intersect_samples,]
sampleSheet_w17 <- sampleSheet_w17[sampleSheet_w17$mouse.id %in% intersect_samples,]
sampleSheet_w24 <- sampleSheet_w24[sampleSheet_w24$mouse.id %in% intersect_samples,]










### Get raw counts for each week. 
otu_count <- otu_count %>%
  remove_rownames() %>%
  column_to_rownames("Taxon") %>%
  `rownames<-`(paste0('OTU-', gsub('usearchOTU','',rownames(.))))


otu_w6_count  <- otu_count %>% select(sampleSheet_w6$original.id)  %>% `colnames<-`(sampleSheet_w6$mouse.id[match(colnames(.),  sampleSheet_w6$original.id)])  %>% t(.)
otu_w17_count <- otu_count %>% select(sampleSheet_w17$original.id) %>% `colnames<-`(sampleSheet_w17$mouse.id[match(colnames(.), sampleSheet_w17$original.id)]) %>% t(.)
otu_w24_count <- otu_count %>% select(sampleSheet_w24$original.id) %>% `colnames<-`(sampleSheet_w24$mouse.id[match(colnames(.), sampleSheet_w24$original.id)]) %>% t(.)



keep_w6  <- unname(apply(otu_w6_count,  2, function(x) mean(x > 0) > 0.25))
keep_w17 <- unname(apply(otu_w17_count, 2, function(x) mean(x > 0) > 0.25))
keep_w24 <- unname(apply(otu_w24_count, 2, function(x) mean(x > 0) > 0.25))



otu_w6_count  <- otu_w6_count[,keep_w6]
otu_w17_count <- otu_w17_count[, keep_w17]
otu_w24_count <- otu_w24_count[, keep_w24]


overlap_otu <- intersect(intersect(colnames(otu_w6_count), colnames(otu_w17_count)), colnames(otu_w24_count))


otu_w6_count  <- otu_w6_count[, overlap_otu]
otu_w17_count <- otu_w17_count[, overlap_otu]
otu_w24_count <- otu_w24_count[, overlap_otu]

taxa_table <- taxa_table[taxa_table$OTU %in% colnames(otu_w6_count),]

otu_w6_count  <- otu_w6_count[, taxa_table$OTU]
otu_w17_count <- otu_w17_count[, taxa_table$OTU]
otu_w24_count <- otu_w24_count[, taxa_table$OTU]


stopifnot(rownames(otu_w6_count) == rownames(otu_w17_count))
stopifnot(rownames(otu_w6_count) == rownames(otu_w24_count))
stopifnot(colnames(otu_w6_count) == colnames(otu_w17_count))
stopifnot(colnames(otu_w6_count) == colnames(otu_w24_count))






### DESeq2 VST
form <- formula(~ 1)
dds_w6  <- DESeqDataSetFromMatrix(countData = t(otu_w6_count),  colData  = sampleSheet_w6, design = form) 
dds_w17 <- DESeqDataSetFromMatrix(countData = t(otu_w17_count), colData = sampleSheet_w17, design = form) 
dds_w24 <- DESeqDataSetFromMatrix(countData = t(otu_w24_count), colData = sampleSheet_w24, design = form) 

otu_vst_w6  <- t(as.matrix(assay(varianceStabilizingTransformation(dds_w6))))
otu_vst_w17 <- t(as.matrix(assay(varianceStabilizingTransformation(dds_w17))))
otu_vst_w24 <- t(as.matrix(assay(varianceStabilizingTransformation(dds_w24))))






### Log2fold change
otu_w6_w17_change  <- log2(otu_w6_count + 1) - log2(otu_w17_count + 1)
otu_w6_w24_change  <- log2(otu_w6_count + 1) - log2(otu_w24_count + 1)
otu_w17_w24_change <- log2(otu_w17_count + 1) - log2(otu_w24_count + 1)










### Rank Z transformation function
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()


otu_w6_w17_change_rz  <- apply(otu_w6_w17_change, 2, rankZ)
otu_w6_w24_change_rz  <- apply(otu_w6_w24_change, 2, rankZ)
otu_w17_w24_change_rz <- apply(otu_w17_w24_change, 2, rankZ)
otu_w6_vst_rz  <- apply(otu_vst_w6,  2, rankZ)
otu_w17_vst_rz <- apply(otu_vst_w17, 2, rankZ)
otu_w24_vst_rz <- apply(otu_vst_w24, 2, rankZ)











### Creating covariate matrix
stopifnot(sampleSheet_w17$diet == sampleSheet_w24$diet)
change.covar <- model.matrix(~diet, data = sampleSheet_w17)[,-1, drop = FALSE]
rownames(change.covar) <- sampleSheet_w17$mouse.id


covar.info <- data.frame(sample.column = 'diet',
                         covar.column  = 'dietprotein',
                         display.name  = 'Diet',
                         interactive   = TRUE,
                         primary       = 'diet',
                         lod.peaks     = 'diet')













## Creating annot.phenotype dataframe
annot.phenotype <- data.frame(data.name   = c(colnames(sampleSheet_w17), colnames(otu_w17_count)),
                              short.name  = c('Mouse Id', 'Coat Color', 'Sex', 'Diet', 'Weel', 'Generation', 'Age at Sacrifice', 'Original Count ID', colnames(otu_w17_count)),
                              R.name      = c(colnames(sampleSheet_w17), colnames(otu_w17_count)),
                              description = c('Mouse identifier', 'Coat color of mouse', 'Sex of mouse: Female (F)', 
                                              'Diet of mouse after week 6: cholic or protein', 'Week at stool sample collection', 'Generation of mouse', 'Age at sacrifice', 'Original mouse identifier in counts matrix given by Weinstock lab',
                                              apply(taxa_table, 1, paste, collapse = '-')),
                              units       = NA,
                              category    = c(rep('demographic', ncol(sampleSheet_w17)), rep('OTU Abundance', ncol(otu_w17_count))),
                              R.category  = c(rep('demographic', ncol(sampleSheet_w17)), rep('OTU Abundance', ncol(otu_w17_count))),
                              is.id       = c(TRUE, rep(FALSE, ncol(sampleSheet_w17) - 1 + ncol(otu_w17_count))),
                              is.numeric  = c(rep(FALSE, ncol(sampleSheet_w17)), rep(TRUE, ncol(otu_w17_count))),
                              is.date     = FALSE,
                              is.factor   = c(rep(FALSE, 3), TRUE, rep(FALSE, ncol(sampleSheet_w17) - 4 + ncol(otu_w17_count))),
                              factor.levels = c(rep(NA, 3), paste0(levels(sampleSheet_w17$diet), collapse = ':'), rep(NA, ncol(sampleSheet_w17) - 4 + ncol(otu_w17_count))),
                              is.covar    = c(rep(FALSE, 3), TRUE, rep(FALSE, ncol(sampleSheet_w17) - 4 + ncol(otu_w17_count))),
                              is.pheno    = c(rep(FALSE, ncol(sampleSheet_w17)), rep(TRUE, ncol(otu_w17_count))),
                              is.derived  = FALSE,
                              omit        = c(rep(FALSE, 4), TRUE, rep(FALSE, ncol(sampleSheet_w17) - 5 + ncol(otu_w17_count))),
                              use.covar   = c(rep(NA, ncol(sampleSheet_w17)), rep('diet', ncol(otu_w17_count))))















### Qtl2 viewer format
dataset.otu.w6 <- list(annot.phenotype = as_tibble(annot.phenotype),
                       annot.samples   = as_tibble(sampleSheet_w6),
                       covar.matrix    = NULL,
                       covar.info      = as_tibble(covar.info),
                       data = list(raw = otu_w6_count,
                                   norm = otu_vst_w6,
                                   rz  = otu_w6_vst_rz),
                       datatype = 'phenotype',
                       display.name = 'Pomp 16s Microbiome: OTU Week 6',
                       lod.peaks = list())




dataset.otu.w17 <- list(annot.phenotype = as_tibble(annot.phenotype),
                        annot.samples   = as_tibble(sampleSheet_w17),
                        covar.matrix    = change.covar,
                        covar.info      = as_tibble(covar.info),
                        data = list(raw = otu_w17_count,
                                    norm = otu_vst_w17,
                                    rz   = otu_w17_vst_rz),
                        datatype = 'phenotype',
                        display.name = 'Pomp 16s Microbiome: OTU Week 17',
                        lod.peaks = list())


dataset.otu.w24 <- list(annot.phenotype = as_tibble(annot.phenotype),
                        annot.samples   = as_tibble(sampleSheet_w24),
                        covar.matrix    = change.covar,
                        covar.info      = as_tibble(covar.info),
                        data = list(raw = otu_w24_count,
                                    norm = otu_vst_w24,
                                    rz   = otu_w24_vst_rz),
                        datatype = 'phenotype',
                        display.name = 'Pomp 16s Microbiome: OTU Week 24',
                        lod.peaks = list())




dataset.otu.w6.w17.change <- list(annot.phenotype = as_tibble(annot.phenotype),
                                  annot.samples   = as_tibble(sampleSheet_w17),
                                  covar.matrix    = change.covar,
                                  covar.info      = as_tibble(covar.info),
                                  data = list(raw.w17  = otu_w6_count,
                                              raw.w24  = otu_w17_count,
                                              log2fold = otu_w6_w17_change,
                                              rz       = otu_w6_w17_change_rz),
                                  datatype = 'phenotype',
                                  display.name = 'Pomp 16s Microbiome: OTU Week 6 - Week 17 Change',
                                  lod.peaks = list())



dataset.otu.w6.w24.change <- list(annot.phenotype = as_tibble(annot.phenotype),
                                  annot.samples   = as_tibble(sampleSheet_w24),
                                  covar.matrix    = change.covar,
                                  covar.info      = as_tibble(covar.info),
                                  data = list(raw.w17  = otu_w6_count,
                                              raw.w24  = otu_w24_count,
                                              log2fold = otu_w6_w24_change,
                                              rz       = otu_w6_w24_change_rz),
                                  datatype = 'phenotype',
                                  display.name = 'Pomp 16s Microbiome: OTU Week 6 - Week 24 Change',
                                  lod.peaks = list())


dataset.otu.w17.w24.change <- list(annot.phenotype = as_tibble(annot.phenotype),
                                   annot.samples   = as_tibble(sampleSheet_w17),
                                   covar.matrix    = change.covar,
                                   covar.info      = as_tibble(covar.info),
                                   data = list(raw.w17  = otu_w17_count,
                                               raw.w24  = otu_w24_count,
                                               log2fold = otu_w17_w24_change,
                                               rz       = otu_w17_w24_change_rz),
                                   datatype = 'phenotype',
                                   display.name = 'Pomp 16s Microbiome: OTU Week 17 - Week 24 Change',
                                   lod.peaks = list())





rm(list = ls()[!ls() %in% c(grep('dataset[.]', ls(), value = TRUE),
                            'genoprobs', 'map', 'markers', 'K')])











### Save
save.image(file = '~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/16S/log2fold/weinstock_pomp_16s_microbiome_qtl_viewer_v2.RData')
