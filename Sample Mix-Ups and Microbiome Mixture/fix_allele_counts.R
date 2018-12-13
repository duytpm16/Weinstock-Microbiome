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
         if(i != 193){
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
      }
      return(results)
}



sample_results_wk6  <- readRDS('Weinstock Microbiome/Mixups/sample_results_wk6_all.rds')
sample_results_wk17 <- readRDS('Weinstock Microbiome/Mixups/sample_results_wk17_all.rds')
sample_results_wk24 <- readRDS('Weinstock Microbiome/Mixups/sample_results_wk24_all.rds')

paired_results_wk6  <- readRDS('Weinstock Microbiome/Mixups/paired_results_wk6_all.rds')
paired_results_wk17 <- readRDS('Weinstock Microbiome/Mixups/paired_results_wk17_all.rds')
paired_results_wk24 <- readRDS('Weinstock Microbiome/Mixups/paired_results_wk24_all.rds')



sample_results_wk6_v2  <- move_AB_to_BB(sample_results_wk6)
sample_results_wk17_v2 <- move_AB_to_BB(sample_results_wk17)
sample_results_wk24_v2 <- move_AB_to_BB(sample_results_wk24)


paired_results_wk6_v2 <- paired_move_AB_to_BB(paired_results_wk6)
paired_results_wk17_v2 <- paired_move_AB_to_BB(paired_results_wk17)
paired_results_wk24_v2 <- paired_move_AB_to_BB(paired_results_wk24)




saveRDS(sample_results_wk6_v2, file = 'Weinstock Microbiome/Mixups/sample_results_wk6_all_v2.rds')
saveRDS(sample_results_wk17_v2, file = 'Weinstock Microbiome/Mixups/sample_results_wk17_all_v2.rds')
saveRDS(sample_results_wk24_v2, file = 'Weinstock Microbiome/Mixups/sample_results_wk24_all_v2.rds')
saveRDS(paired_results_wk6_v2, file = 'Weinstock Microbiome/Mixups/paired_results_wk6_all_v2.rds')
saveRDS(paired_results_wk17_v2, file = 'Weinstock Microbiome/Mixups/paired_results_wk17_all_v2.rds')
paired_results_wk24_v2 <- paired_results_wk24_v2[-193]
saveRDS(paired_results_wk24_v2, file = 'Weinstock Microbiome/Mixups/paired_results_wk24_all_v2.rds')

