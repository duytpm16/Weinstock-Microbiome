#####################################################################################################################
#
#   This script is used to edit the POMP genoprobs after downloading from dodb.jax.org
#
#
#
#
#   Input:
#     1.) MegaMuga SNPs Annotation
#     2.) Qtl2 genotype probability array from dodb.jax.org
#
#
#   Output:
#     1.) New .RData file containing new genoprobs, K, map, and markers
#
#
#
#
#   Author: Duy Pham
#   Date:   April 12th, 2019
#   E-mail: duy.pham@jax.org
#####################################################################################################################

### Options and Libraries
options(stringsAsFactors = FALSE)
library(dplyr)
library(qtl2)
library(qtl2convert)









### Load data
load("/Users/phamd/Downloads/snps.megamuga.Rdata")
genoprobs <- readRDS("~/Desktop/Weinstock Pomp Microbiome/Pomp-194_Pomp_DO/Pomp_194_Pomp_DO__genoprobs_8state_MegaMUGA.Rdata")


















### Editing names of markers data frame and genoprobabilities. There are two duplicates: DO-420 and DO-595. Adding a .1 to the sample that has a lower call rate
snps <- snps %>%
             mutate(chr = gsub('chr', '', .$chr),
                    chr = factor(chr, levels = c(1:19,'X')),
                    bp  = pos,
                    pos = bp / 1e6) %>%
             filter(chr %in% c(1:19,'X')) %>%
             arrange(chr, cM, pos) %>%
             dplyr::select(-A1, -A2, -seq.A, -seq.B)





genoprobs <- probs_qtl2_to_doqtl(probs = genoprobs)
genoprobs <- genoprobs[grep('DP-DO2', dimnames(genoprobs)[[1]]),,]
dimnames(genoprobs)[[1]] <- gsub('UNC-Pomp_Mouse_12dec2012_DP-', '', dimnames(genoprobs)[[1]])
dimnames(genoprobs)[[1]] <- gsub("\\_.*", '', dimnames(genoprobs)[[1]])
dimnames(genoprobs)[[1]] <- gsub("-F", '', dimnames(genoprobs)[[1]])
dimnames(genoprobs)[[1]][duplicated(dimnames(genoprobs)[[1]], fromLast = TRUE)] <- paste0(dimnames(genoprobs)[[1]][duplicated(dimnames(genoprobs)[[1]], fromLast = TRUE)], '.1')
genoprobs <- genoprobs[order(dimnames(genoprobs)[[1]]),,]















### Qtl2 data
genoprobs <- probs_doqtl_to_qtl2(probs = genoprobs,
                                 map   = snps,
                                 chr_column = 'chr',
                                 pos_column = 'pos',
                                 marker_column = 'marker')



K <- calc_kinship(probs = genoprobs,
                  type  = 'loco',
                  cores = 0)



markers <- snps %>%
                mutate(chr = as.character(chr))

map     <- map_df_to_list(map = markers,
                          chr_column =  'chr',
                          pos_column =  'pos',
                          marker_column = 'marker')

rm(snps)
