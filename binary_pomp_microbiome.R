##################################################################################################################
#   Try to process qtl scans with microbiome data as binary. 0 - not present          1 - present.
#
#   Input:
#       1.) .RData file generated from gather_
#
#
#   Output:
#       1.) New .RData file without rZ data and microbiome data as binary
#
#
#   Author: Duy Pham
#   E-mail: duy.pham@jax.org
#   Date: October 23, 2018
##################################################################################################################
load("~/Desktop/Weinstock Microbiome/weinstock_16s_microbiome_qtl2_input.RData")



### Save data as binary matrix
for(i in ls()[grep('_w6', ls())]){
    
    if(!(grepl('samples', i) || grepl('rZ', i) || grepl('covar',i))){
       temp <- get(i)
       temp[temp > 0 ] <- 1
       assign(i, temp)
    }
}


for(i in ls()[grep('_w17', ls())]){
  
  if(!(grepl('samples', i) || grepl('rZ', i) || grepl('covar',i))){
    temp <- get(i)
    temp[temp > 0 ] <- 1
    assign(i, temp)
  }
}



for(i in ls()[grep('_w24', ls())]){
  
  if(!(grepl('samples', i) || grepl('rZ', i) || grepl('covar',i))){
    temp <- get(i)
    temp[temp > 0 ] <- 1
    assign(i, temp)
  }
}

rm(list = c(ls()[grep('rZ',ls())]))
rm(temp, i)



save.image('weinstock_16s_microbiome_binary_qtl2_input.RData')
