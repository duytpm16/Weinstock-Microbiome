move_AB_to_BB <-  function(results){
  
  for(i in 1:length(results)){
    
    for(j in 1:dim(results[[i]])[1]){
      
      AB <- results[[i]][j,,][2,]
      results[[i]][j,,][3,] <- AB
      
      results[[i]][j,,][2,] <- c(0,0)
    }
  }
  return(results)
}



paired_move_AB_to_BB <-  function(results){
  
  for(i in 1:length(results)){
      for(k in 1:2){
        for(j in 1:dim(results[[i]])[1]){
          
          AA.BB <- results[[i]][,,,k][j,,][1,2]
          BB.AA <- results[[i]][,,,k][j,,][2,1]
          BB.BB <- results[[i]][,,,k][j,,][2,2]
          
          results[[i]][,,,k][j,,][1,3] <- AA.BB
          results[[i]][,,,k][j,,][3,1] <- BB.AA
          results[[i]][,,,k][j,,][3,3] <- BB.BB
          
          
          results[[i]][,,,k][j,,][1,2] <- 0
          results[[i]][,,,k][j,,][2,1] <- 0
          results[[i]][,,,k][j,,][2,2] <- 0
      }
    }
  }
  return(results)
}



### Read in data
sample_results_wk6  <- readRDS('sample_results_all_wk_6.rds')
sample_results_wk17 <- readRDS('sample_results_all_wk_17.rds')
sample_results_wk24 <- readRDS('sample_results_all_wk_24.rds')

paired_results_wk6  <- readRDS('paired_results_all_wk_6.rds')
paired_results_wk17 <- readRDS('paired_results_all_wk_17.rds')
paired_results_wk24 <- readRDS('paired_results_all_wk_24.rds')




### Removing these samples because there is 0 counts
rm_wk6_samples  <- which(names(sample_results_wk6) %in% c('DPDP.DO2.573.F', 'DPDP.DO2.574.F', 'DPDP.DO2.575.F', 'DPDP.DO2.811.F'))
rm_wk17_samples <- which(names(sample_results_wk17) %in% c('DPDP.DO2.742.F', 'DPDP.DO2.746.F', 'DPDP.DO2.744.F'))
rm_wk24_samples <- which(names(sample_results_wk24) %in% c('DPDP.DO2.500.F', 'DPDP.DO2.746.F', 'DPDP.DO2.811.F'))
rm_wk6_paired_samples  <- which(names(paired_results_wk6) %in% c('DPDP.DO2.513.F', 'DPDP.DO2.533.F', 'DPDP.DO2.548.F', 
                                                                 'DPDP.DO2.573.F','DPDP.DO2.574.F','DPDP.DO2.575.F','DPDP.DO2.591.F',
                                                                 'DPDP.DO2.593.F','DPDP.DO2.675.F','DPDP.DO2.703.F','DPDP.DO2.711.F',
                                                                 'DPDP.DO2.837.F','DPDP.DO2.838.F','DPDP.DO2.811.F'))
rm_wk17_paired_samples <- which(names(paired_results_wk17) %in% c('DPDP.DO2.513.F', 'DPDP.DO2.533.F', 'DPDP.DO2.548.F', 
                                                                  'DPDP.DO2.591.F', 'DPDP.DO2.593.F','DPDP.DO2.675.F',
                                                                  'DPDP.DO2.703.F','DPDP.DO2.711.F','DPDP.DO2.739.F','DPDP.DO2.742.F',
                                                                  'DPDP.DO2.744.F', 'DPDP.DO2.746.F', 'DPDP.DO2.837.F','DPDP.DO2.838.F'))
rm_wk24_paired_samples <- which(names(paired_results_wk24) %in% c('DPDP.DO2.500.F', 'DPDP.DO2.513.F', 'DPDP.DO2.533.F', 'DPDP.DO2.548.F',
                                                                  'DPDP.DO2.591.F', 'DPDP.DO2.593.F','DPDP.DO2.675.F',
                                                                  'DPDP.DO2.703.F','DPDP.DO2.711.F','DPDP.DO2.746.F', 
                                                                  'DPDP.DO2.811.F','DPDP.DO2.837.F','DPDP.DO2.838.F'))


sample_results_wk6  <- sample_results_wk6[-rm_wk6_samples]
sample_results_wk17 <- sample_results_wk17[-rm_wk17_samples]
sample_results_wk24 <- sample_results_wk24[-rm_wk24_samples]

paired_results_wk6  <- paired_results_wk6[-rm_wk6_paired_samples]
paired_results_wk17 <- paired_results_wk17[-rm_wk17_paired_samples]
paired_results_wk24 <- paired_results_wk24[-rm_wk24_paired_samples]






### Moving AB counts to BB counts
sample_results_wk6_v2  <- move_AB_to_BB(sample_results_wk6)
sample_results_wk17_v2 <- move_AB_to_BB(sample_results_wk17)
sample_results_wk24_v2 <- move_AB_to_BB(sample_results_wk24)

paired_results_wk6_v2  <- paired_move_AB_to_BB(paired_results_wk6)
paired_results_wk17_v2 <- paired_move_AB_to_BB(paired_results_wk17)
paired_results_wk24_v2 <- paired_move_AB_to_BB(paired_results_wk24)








### Save file
saveRDS(sample_results_wk6_v2, file = 'Weinstock Microbiome/Mixups/sample_results_wk6_all_v2.rds')
saveRDS(sample_results_wk17_v2, file = 'Weinstock Microbiome/Mixups/sample_results_wk17_all_v2.rds')
saveRDS(sample_results_wk24_v2, file = 'Weinstock Microbiome/Mixups/sample_results_wk24_all_v2.rds')
saveRDS(paired_results_wk6_v2, file = 'Weinstock Microbiome/Mixups/paired_results_wk6_all_v2.rds')
saveRDS(paired_results_wk17_v2, file = 'Weinstock Microbiome/Mixups/paired_results_wk17_all_v2.rds')
saveRDS(paired_results_wk24_v2, file = 'Weinstock Microbiome/Mixups/paired_results_wk24_all_v2.rds')
