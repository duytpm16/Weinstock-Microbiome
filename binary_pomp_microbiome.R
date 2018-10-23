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
