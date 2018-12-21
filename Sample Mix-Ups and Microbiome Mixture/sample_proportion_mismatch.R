samp_wk6  <- readRDS('Weinstock Microbiome/Mixups/sample_results_wk6_all_v2.rds')
samp_wk17 <- readRDS('Weinstock Microbiome/Mixups/sample_results_wk17_all_v2.rds')
samp_wk24 <- readRDS('Weinstock Microbiome/Mixups/sample_results_wk24_all_v2.rds')



# calculate proportion of mismatches at homozygous loci
f <- function(a, rn) {
  x <- apply(a, 1, function(b) (b[1,2] + b[3,1]) / sum(b[1,] + b[3,]))
  x[rn]
}



samp_prop_mis_wk6  <- t(sapply(samp_wk6,  f, rownames(samp_wk6[[1]])))
samp_prop_mis_wk17 <- t(sapply(samp_wk17, f, rownames(samp_wk17[[1]])))
samp_prop_mis_wk24 <- t(sapply(samp_wk24, f, rownames(samp_wk24[[1]])))





saveRDS(samp_prop_mis_wk6,  file="Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk6.rds")
saveRDS(samp_prop_mis_wk17, file="Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk17.rds")
saveRDS(samp_prop_mis_wk24, file="Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk24.rds")