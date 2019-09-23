### Options and Libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(qtl2)









### Load data
load("~/Desktop/Weinstock_Pomp_Microbiome/Genotypes/Original Genotypes/weinstock_pomp_genoprobs_map_K_qtl2.Rdata")
metadata <- read.csv("~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/Mouse Info/DO1DO2micecombined_Final.csv")
sampleSheet <- read.xlsx(xlsxFile = "~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/16S/Raw Counts/taxonomic_abundance_filtered_unnormalized_20181004.xlsx")
taxa_table  <- read.xlsx(xlsxFile = "~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/16S/Raw Counts/taxonomic_abundance_filtered_unnormalized_20181004.xlsx", sheet = 'OTU Taxonomy')
genus_count <- read.xlsx(xlsxFile = "~/Desktop/Weinstock_Pomp_Microbiome/Phenotypes/16S/Raw Counts/taxonomic_abundance_filtered_unnormalized_20181004.xlsx", sheet = 'genus_abundance')

bad_samples <- c('DO2.444','DO2.492','DO2.507','DO2.552','DO2.541','DO2.440',
                 'DO2.484','DO2.833','DO2.835','DO2.836','DO2.840','DO2.841',
                 'DO2.850','DO2.411','DO2.416','DO2.472','DO2.504','DO2.539',
                 'DO2.549','DO2.732','DO2.801','DO2.825','DO2.846')














### Editing sample sheet
sampleSheet <- sampleSheet  %>% 
                  mutate(MouseID     = paste0('DO2.', Mouse), 
                         sex         = 'F', 
                         generation  = '11',
                         original.id = paste0('16s-', SampleName)) %>%
                  full_join(y = metadata[,c('MouseID', 'CtClr','AgeSac')], by = 'MouseID', copy = TRUE) %>%
                  filter(grepl('DO2', MouseID)) %>%
                  dplyr::rename(mouse.id = MouseID, coat.color = CtClr, age.sac = AgeSac, diet = Diet, week = Week) %>%
                  mutate(coat.color = case_when(coat.color == 1 ~ 'white', coat.color == 2 ~ 'agouti', coat.color == 3 ~ 'black', coat.color == 4 ~ 'other')) %>%
                  distinct() %>% arrange(diet, week) %>%
                  select(mouse.id, sex, diet, week, generation, coat.color, age.sac, original.id)













### Create sample sheet for each week
sampleSheet_w6  <- sampleSheet[sampleSheet$week == 6,]  %>% mutate(sex = factor(sex), diet = factor(diet)) %>% `rownames<-`(.$mouse.id)
sampleSheet_w17 <- sampleSheet[sampleSheet$week == 17,] %>% mutate(sex = factor(sex), diet = factor(diet)) %>% `rownames<-`(.$mouse.id)
sampleSheet_w24 <- sampleSheet[sampleSheet$week == 24,] %>% mutate(sex = factor(sex), diet = factor(diet)) %>% `rownames<-`(.$mouse.id)


samples <- list(w6 = sampleSheet_w6$mouse.id, w17 = sampleSheet_w17$mouse.id, w24 = sampleSheet_w24$mouse.id, gp = dimnames(genoprobs[[1]])[[1]])
samples <- Reduce(intersect, samples)
samples <- samples[!samples %in% bad_samples]


sampleSheet_w6  <- sampleSheet_w6[sampleSheet_w6$mouse.id %in% samples,]
sampleSheet_w17 <- sampleSheet_w17[sampleSheet_w17$mouse.id %in% samples,]
sampleSheet_w24 <- sampleSheet_w24[sampleSheet_w24$mouse.id %in% samples,]

















### Editing taxa table
taxa_table <- apply(taxa_table, 2, FUN = function(x) gsub('"', "",  x, fixed = TRUE))
taxa_table <- apply(taxa_table, 2, FUN = function(x) gsub('_', " ", x, fixed = TRUE))
taxa_table <- taxa_table %>%
                as.data.frame() %>%
                select(genus, family, order, class, phylum, domain) %>%
                arrange(genus) %>%
                distinct()














### Get raw counts for each week. 
genus_count <- genus_count %>%
                  mutate(Taxon = gsub('"', '', Taxon, fixed = TRUE)) %>%
                  mutate(Taxon = gsub('_', ' ', Taxon, fixed = TRUE)) %>%
                  column_to_rownames("Taxon")


genus_w6_count  <- genus_count %>% select(sampleSheet_w6$original.id)  %>% `colnames<-`(sampleSheet_w6$mouse.id[match(colnames(.),  sampleSheet_w6$original.id)])  %>% t(.)
genus_w6_count  <- genus_w6_count[sampleSheet_w6$mouse.id,]

genus_w17_count <- genus_count %>% select(sampleSheet_w17$original.id) %>% `colnames<-`(sampleSheet_w17$mouse.id[match(colnames(.), sampleSheet_w17$original.id)]) %>% t(.)
genus_w17_count <- genus_w17_count[sampleSheet_w17$mouse.id,]

genus_w24_count <- genus_count %>% select(sampleSheet_w24$original.id) %>% `colnames<-`(sampleSheet_w24$mouse.id[match(colnames(.), sampleSheet_w24$original.id)]) %>% t(.)
genus_w24_count <- genus_w24_count[sampleSheet_w24$mouse.id,]












### Reorder counts dataframe
overlap_genus <- intersect(intersect(colnames(genus_w6_count), colnames(genus_w17_count)), colnames(genus_w24_count))
taxa_table  <- taxa_table[taxa_table$genus %in% overlap_genus,]

genus_w6_count  <- genus_w6_count[,  taxa_table$genus]
genus_w17_count <- genus_w17_count[, taxa_table$genus]
genus_w24_count <- genus_w24_count[, taxa_table$genus]












### DESeq2 VST
form <- formula(~ 1)
dds_wk6 <- DESeqDataSetFromMatrix(countData = t(genus_w6_count), colData  = sampleSheet_w6, design = form) 
med_wk6 <- t(counts(estimateSizeFactors(dds_wk6, type = 'poscounts'), normalized = TRUE))
vst_wk6 <- t(as.matrix(assay(varianceStabilizingTransformation(dds_wk6))))


form <- formula(~ diet)
dds_wk17 <- DESeqDataSetFromMatrix(countData = t(genus_w17_count), colData  = sampleSheet_w17, design = form) 
med_wk17 <- t(counts(estimateSizeFactors(dds_wk17, type = 'poscounts'), normalized = TRUE))
vst_wk17 <- t(as.matrix(assay(varianceStabilizingTransformation(dds_wk17))))


form <- formula(~ diet)
dds_wk24 <- DESeqDataSetFromMatrix(countData = t(genus_w24_count), colData  = sampleSheet_w24, design = form) 
med_wk24 <- t(counts(estimateSizeFactors(dds_wk24, type = 'poscounts'), normalized = TRUE))
vst_wk24 <- t(as.matrix(assay(varianceStabilizingTransformation(dds_wk24))))









### Rank Z transformation function
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()


rz_wk6  <- apply(vst_wk6,  2, rankZ)
rz_wk17 <- apply(vst_wk17, 2, rankZ)
rz_wk24 <- apply(vst_wk24, 2, rankZ)












### Creating covariate matrix
covar.w6  <- NULL
covar.w17 <- model.matrix(~ diet, data = sampleSheet_w17)[, -1, drop = FALSE]
colnames(covar.w17) <- 'diet'
covar.w24 <- model.matrix(~ diet, data = sampleSheet_w24)[, -1, drop = FALSE]
colnames(covar.w24) <- 'diet'









### Creating covar.info dataframe for qtl2 viewer
covar.info <- data.frame(sample.column = 'diet',
                         covar.column  = 'diet',
                         display.name  = 'Diet',
                         interactive   = TRUE,
                         primary       = 'diet',
                         lod.peaks     = 'diet_int')













## Creating annot.phenotype dataframe
annot.phenotype <- data.frame(data.name   = c(colnames(sampleSheet), taxa_table$genus),
                              short.name  = c('Mouse Id', 'Sex', 'Diet', 'Week', 'Generation', 'Coat Color', 'Age at Sacrifice', 'Original Count ID', taxa_table$genus),
                              R.name      = c(colnames(sampleSheet), taxa_table$genus),
                              description = c('Mouse identifier', 'Sex of mouse: Female (F) - Male (M)', 'Diet of mouse', 
                                              'Week at stool sample collection', 'Generation of mouse', 'Coat color of mouse', 'Age at sacrifice', 'Original mouse identifier in counts matrix given by Weinstock lab',
                                              apply(taxa_table, 1, paste, collapse = ' - ')),
                              units       = NA,
                              category    = c(rep('Demographic', ncol(sampleSheet)), rep('genus Abundance', length(taxa_table$genus))),
                              R.category  = c(rep('Demographic', ncol(sampleSheet)), rep('genus Abundance', length(taxa_table$genus))),
                              is.id       = c(TRUE, rep(FALSE, ncol(sampleSheet) - 1 + length(taxa_table$genus))),
                              is.numeric  = c(rep(FALSE, ncol(sampleSheet)), rep(TRUE, length(taxa_table$genus))),
                              is.date     = FALSE,
                              is.factor   = c(FALSE, FALSE, TRUE, rep(FALSE, ncol(sampleSheet) - 3 + length(taxa_table$genus))),
                              factor.levels = c(rep(NA, 2), paste0(levels(sampleSheet_w17$diet), collapse = ':'), rep(NA, ncol(sampleSheet) - 3 + length(taxa_table$genus))),
                              is.covar    = c(rep(FALSE, 2), TRUE, rep(FALSE, ncol(sampleSheet_w17) - 3 + length(taxa_table$genus))),
                              is.pheno    = c(rep(FALSE, ncol(sampleSheet)), rep(TRUE, length(taxa_table$genus))),
                              is.derived  = FALSE,
                              omit        = FALSE,
                              use.covar   = c(rep(NA, ncol(sampleSheet)), rep('diet', length(taxa_table$genus))))















### Qtl2 viewer format
dataset.genus.w6 <- list(annot.phenotype = as_tibble(annot.phenotype),
                       annot.samples   = as_tibble(sampleSheet_w6),
                       covar.matrix    = covar.w6,
                       covar.info      = NULL,
                       data = list(raw = genus_w6_count,
                                   vst = vst_wk6,
                                   rz  = rz_wk6),
                       taxa = taxa_table,
                       datatype = 'phenotype',
                       display.name = 'Genus - Week 6',
                       lod.peaks = list())




dataset.genus.w17 <- list(annot.phenotype = as_tibble(annot.phenotype),
                        annot.samples   = as_tibble(sampleSheet_w17),
                        covar.matrix    = covar.w17,
                        covar.info      = covar.info,
                        data = list(raw = genus_w17_count,
                                    vst = vst_wk17,
                                    rz  = rz_wk17),
                        taxa = taxa_table,
                        datatype     = 'phenotype',
                        display.name = 'Genus - Week 17',
                        lod.peaks    = list())





dataset.genus.w24 <- list(annot.phenotype = as_tibble(annot.phenotype),
                        annot.samples   = as_tibble(sampleSheet_w24),
                        covar.matrix    = covar.w24,
                        covar.info      = covar.info,
                        data = list(raw = genus_w24_count,
                                    vst = vst_wk24,
                                    rz  = rz_wk24),
                        taxa = taxa_table,
                        datatype     = 'phenotype',
                        display.name = 'Genus - Week 24',
                        lod.peaks    = list())










### Save
rm(list = ls()[!ls() %in% c(grep('dataset[.]', ls(), value = TRUE))])
load('~/Desktop/weinstock_pomp_16s_microbiome_qtl_viewer_v4.RData')
save.image(file = '~/Desktop/weinstock_pomp_16s_microbiome_qtl_viewer_v4.RData')
