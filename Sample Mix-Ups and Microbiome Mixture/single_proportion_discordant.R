single_results_wk6_all_v2  <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk6_all_v2.rds")
single_results_wk17_all_v2 <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk17_all_v2.rds")
single_results_wk24_all_v2 <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk24_all_v2.rds")





# calculate proportion of mismatches at homozygous loci
f <- function(a, rn ) {
  x <- apply(a, 1, function(b) (b[1,2] + b[3,1]) / sum(b[1,] + b[3,]))
  x[rn] }





samp_p_w6 <- t(sapply(single_results_wk6_all_v2, f, rownames(single_results_wk6_all_v2[[1]])))
samp_p_w17 <- t(sapply(single_results_wk17_all_v2, f, rownames(single_results_wk17_all_v2[[1]])))
samp_p_w24 <- t(sapply(single_results_wk24_all_v2, f, rownames(single_results_wk24_all_v2[[1]])))
saveRDS(samp_p, file="single_results.rds")