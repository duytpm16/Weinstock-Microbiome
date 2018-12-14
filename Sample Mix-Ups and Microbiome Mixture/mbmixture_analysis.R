options(stringsAsFactors = FALSE)
library(dplyr)
library(mbmixture)
library(parallel)



args  <- commandArgs(trailingOnly = TRUE)
pair  <- readRDS(args[1])
week  <- args[2]
cores <- as.numeric(args[3])


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




saveRDS(result, file=paste0("mixture_results_wk_",week,".rds"))

