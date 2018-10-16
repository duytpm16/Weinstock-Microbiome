options(stringsAsFactors = FALSE)

library(qtl2)


### Load in qtl scans and map
genus_qtl <- readRDS('weinstock_16s_microbiome_genus_rZ_qtl_lod.rds')
family_qtl <- readRDS('weinstock_16s_microbiome_family_rZ_qtl_lod.rds')
class_qtl <- readRDS('weinstock_16s_microbiome_class_rZ_qtl_lod.rds')
order_qtl <- readRDS('weinstock_16s_microbiome_order_rZ_qtl_lod.rds')
phylum_qtl <- readRDS('weinstock_16s_microbiome_phylum_rZ_qtl_lod.rds')


load("~/Desktop/weinstock_16s_microbiome_qtl2_input.RData")





### Function to create lod peaks table with BLUP scans
findpeaks_table <- function(expr, covar, genoprobs, kinship, qtl_output, map, threshold){
  
  
  
      # QTL2 find_peaks function
      lod.peaks <- find_peaks(qtl_output, map = map, threshold = threshold)
      
      
  
  
      
      # Changing colnames of find_peaks output
      marker.id <- paste0(as.character(lod.peaks$chr), '_', round(lod.peaks$pos * 1000000))
      annot.id <- lod.peaks[,'lodcolumn']
      lod.peaks <- cbind(annot.id, marker.id, lod.peaks[,c('lod','chr','pos')])
      colnames(lod.peaks) <- c('annot.id','marker.id','lod','qtl.chr','qtl.pos')
      lod.peaks$qtl.chr <- as.character(lod.peaks$qtl.chr)
      
      
  
  
      
      # BLUP scan
      lod.peaks = cbind(lod.peaks, matrix(0, nrow = nrow(lod.peaks), ncol = 8, 
                                          dimnames = list(NULL, LETTERS[1:8])))
      for(i in 1:nrow(lod.peaks)) {
        
          chr  = lod.peaks$qtl.chr[i]
          mkr  = lod.peaks$marker.id[i]
          annot = lod.peaks$annot.id[i]
        
          # BLUP scan at QTL.
          gp = genoprobs[,chr]
          gp[[1]] = gp[[1]][,,mkr,drop=FALSE]
          blup = scan1blup(genoprobs = gp, pheno = expr[, annot, drop = FALSE],
                           kinship = K[[chr]], addcovar = covar, cores = 0)
        
          lod.peaks[i,6:13] = blup[1,1:8]
      }
      
      
      
      
      return(lod.peaks)
}










### Get all taxa qtl scans and add a lod peaks dataframe where peaks are > 6 for each week.
for(i in ls()[grep('qtl',ls())]){
    temp <- get(i)

    temp[["peaks_w6"]] <- findpeaks_table(expr = get(paste0('rZ_', gsub('_qtl', '_w6',  i))), 
                                          qtl_output = temp$w6, covar = covar_w6, genoprobs = genoprobs, 
                                          kinship = K, map = map, threshold = 6)
    
    temp[["peaks_w17"]] <- findpeaks_table(expr = get(paste0('rZ_', gsub('_qtl', '_w17',  i))), 
                                           qtl_output = temp$w17, covar = covar_w17, genoprobs = genoprobs, 
                                           kinship = K, map = map, threshold = 6)
    
    temp[["peaks_w24"]] <- findpeaks_table(expr = get(paste0('rZ_', gsub('_qtl', '_w24',  i))), 
                                           qtl_output = temp$w24, covar = covar_w24, genoprobs = genoprobs, 
                                           kinship = K, map = map, threshold = 6)
    
    assign(i, temp)
}








saveRDS(genus_qtl, 'weinstock_16s_microbiome_genus_rZ_qtl.rds')
saveRDS(family_qtl, 'weinstock_16s_microbiome_family_rZ_qtl.rds')
saveRDS(class_qtl, 'weinstock_16s_microbiome_class_rZ_qtl.rds')
saveRDS(order_qtl, 'weinstock_16s_microbiome_order_rZ_qtl.rds')
saveRDS(phylum_qtl, 'weinstock_16s_microbiome_phylum_rZ_qtl.rds')


