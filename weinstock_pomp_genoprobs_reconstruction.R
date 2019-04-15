### Options and Libraries
options(stringsAsFactors = FALSE)
library(dplyr)
library(qtl2)
library(qtl2convert)









### Load Data
pheno <- read_csv('~/Desktop/Weinstock Pomp Microbiome/Phenotypes/DO1DO2micecombined_Final.csv')
geno  <- read_csv('~/Desktop/Weinstock Pomp Microbiome/Genotypes/Pomp-194_Pomp_DO-MegaMUGA/Pomp-194_Pomp_DO-MegaMUGA_geno.csv')








### Subset and create covar
pheno <- pheno %>%
               filter(DO == 2) %>%
               mutate(MouseID = paste0('DP-', gsub('.','-',MouseID, fixed = TRUE), '-F'))

geno  <- geno[,c('marker', grep('DO2', colnames(geno), value = TRUE))]

covar <- data.frame(id  = colnames(geno),
                    Sex = rep('F', ncol(geno)),
                    gen = rep('11', ncol(geno)))








### Write new geno and covar .csv files
write.csv(x = geno, file = '~/Desktop/Weinstock Pomp Microbiome/Genotypes/Pomp-194_Pomp_DO-MegaMUGA/Pomp-194_Pomp_DO-MegaMUGA_geno.csv', row.names = FALSE)
write.csv(x = covar, file = '~/Desktop/Weinstock Pomp Microbiome/Genotypes/Pomp-194_Pomp_DO-MegaMUGA/_covar.csv', row.names = FALSE)








### Read in cross
pomp_cross <- read_cross2(file = '~/Desktop/Weinstock Pomp Microbiome/Genotypes/Pomp-194_Pomp_DO-MegaMUGA/Pomp-194_Pomp_DO-MegaMUGA.json')








### Replacing the geno, pmap, and gmap in pomp_cross with just the autosome and X chromosome
geno <- list()
pmap <- list()
gmap <- list()

for(i in c(1:19,'X')){
  
  geno[[i]] <- pomp_cross$geno[[i]] 
  pmap[[i]] <- pomp_cross$pmap[[i]]
  gmap[[i]] <- pomp_cross$gmap[[i]]
  
}

pomp_cross$geno <- geno
pomp_cross$pmap <- pmap
pomp_cross$gmap <- gmap








### Calculate genotype probabilities
probs <- calc_genoprob(cross = pomp_cross, cores = 0)



### Convert prob to allele probabilities
genoprobs <- genoprob_to_alleleprob(probs = probs, cores = 0)



### Calculate Kinship
K <- calc_kinship(probs = genoprobs, type = 'loco', cores = 0)












### Get physical map list
map <- pmap



### Get markers dataframe and combine cM
markers   <- map_list_to_df(map)
g_markers <- map_list_to_df(gmap)

markers <- merge(markers, g_markers, by = 'marker', all.x = TRUE, sort = FALSE)
markers <- markers %>%
                   select(marker, chr.x, pos.x, pos.y) %>%
                   dplyr::rename(chr = chr.x, pos = pos.x, cM = pos.y)

rownames(markers) <- markers$marker







save(markers,map,genoprobs, K, file = '~/Desktop/Weinstock Pomp Microbiome/MetaBAT/weinstock_pomp_genoprobs_K_map_markers.RData')