### Load in combined results
results_wk6  <- readRDS("sample_results_week_6_all.rds")
results_wk17 <- readRDS("sample_results_week_17_all.rds")
results_wk24 <- readRDS("sample_results_week_24_all.rds")







### Calculate proportion of mismatches at homozygous loci (Karl's code)
f <- function(a, rn ) {
  x <- apply(a, 1, function(b) (b[1,2] + b[3,1]) / sum(b[1,] + b[3,]))
  x[rn] 
}


w6  <- t(sapply(results_wk6,  f, rownames(results_wk6[[1]])))
w17 <- t(sapply(results_wk17, f, rownames(results_wk17[[1]])))
w24 <- t(sapply(results_wk24, f, rownames(results_wk24[[1]])))


             
             
             
             
             
             
             
### Save
saveRDS(wk6,  file="single_results_prop_mismatch_w6.rds")
saveRDS(wk17, file="single_results_prop_mismatch_w17.rds")
saveRDS(wk24, file="single_results_prop_mismatch_w24.rds")
