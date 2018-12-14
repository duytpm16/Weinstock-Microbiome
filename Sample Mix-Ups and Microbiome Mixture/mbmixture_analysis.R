### Options and Libraries
options(stringsAsFactors = FALSE)
library(dplyr)
library(mbmixture)
library(parallel)





### Read in required data
#     1.) pair:  pair result list
#     2.) week:  which week
#     3.) cores: how many cores
args  <- commandArgs(trailingOnly = TRUE)
pair  <- readRDS(args[1])
week  <- args[2]
cores <- as.numeric(args[3])






### Do MB mixture analysis for each sample in pair
analyze_one <- function(i) {
  res <- matrix(nrow=nrow(pair[[i]]), ncol=5)
  
  for(j in 1:nrow(pair[[i]])) {
    if(names(pair)[i] == rownames(pair[[i]])[j]) next
    res[j,] <- tmp <- mle_pe(pair[[i]][j,,,])
  }
  
  dimnames(res) <- list(rownames(pair[[i]]), names(tmp))
  res
}






result <- mclapply(seq_along(pair), analyze_one, mc.cores=cores)
names(result) <- names(pair)






### Save results
saveRDS(result, file=paste0("mixture_results_wk_",week,".rds"))

