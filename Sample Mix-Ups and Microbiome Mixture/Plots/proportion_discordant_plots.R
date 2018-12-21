library(lineup2)
library(broman)




sample_results_wk6  <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk6_all_v2.rds")
sample_results_wk17 <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk17_all_v2.rds")
sample_results_wk24 <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk24_all_v2.rds")
samp_wk6  <- readRDS("Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk6.rds")
samp_wk17 <- readRDS("Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk17.rds")
samp_wk24 <- readRDS("Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk24.rds")





### Week 6
grayplot(x = get_self(samp_wk6), y = get_best(samp_wk6),
         xlab = "Self-self distance", ylab = "Best distance",
         xlim=c(0, 0.213), ylim=c(0, 0.213), 
         xaxs="i", yaxs="i",
         main = 'Proportion of Mismatch Week 6')
plot_label <- gsub('DO2.','',gsub('.F','', gsub('DPDP.','',rownames(samp_wk6))))
text(get_self(samp_wk6), get_best(samp_wk6) + 0.005, labels = plot_label)




lowcounts_samples <- c('DPDP.DO2.403.F','DPDP.DO2.410.F','DPDP.DO2.413.F','DPDP.DO2.429.F','DPDP.DO2.558.F','DPDP.DO2.567.F','DPDP.DO2.571.F','DPDP.DO2.806.F')
points(get_self(samp_wk6)[lowcounts_samples], get_best(samp_wk6)[lowcounts_samples], pch = 21, bg = 'green')
mixture_samples   <- c('DPDP.DO2.444.F','DPDP.DO2.492.F','DPDP.DO2.552.F')
points(get_self(samp_wk6)[mixture_samples], get_best(samp_wk6)[mixture_samples], pch = 21, bg = 'red')








### Week 17
grayplot(get_self(samp_wk17), get_best(samp_wk17),
         xlab="Self-self distance", ylab="Best distance",
         xlim=c(0, 0.213), ylim=c(0, 0.213), xaxs="i", yaxs="i")
plot_label <- gsub('.F','', gsub('DPDP.','',rownames(samp_wk17)))
text(get_self(samp_wk17), get_best(samp_wk17) + 0.005, labels = plot_label)
 

lowcounts_samples <- c('DPDP.DO2.415.F','DPDP.DO2.741.F','DPDP.DO2.745.F','DPDP.DO2.550.F')
points(get_self(samp_wk17)[lowcounts_samples], get_best(samp_wk17)[lowcounts_samples], pch = 21, bg = 'green')
mixture_samples   <- c('DPDP.DO2.440.F','DPDP.DO2.437.F','DPDP.DO2.484.F','DPDP.DO2.552.F', 'DPDP.DO2.835.F','DPDP.DO2.836.F','DPDP.DO2.850.F')
points(get_self(samp_wk17)[mixture_samples], get_best(samp_wk17)[mixture_samples], pch = 21, bg = 'red')







### Week 24
grayplot(get_self(samp_wk24), get_best(samp_wk24),
         xlab="Self-self distance", ylab="Best distance",
         xlim=c(0, 0.213), ylim=c(0, 0.213), xaxs="i", yaxs="i")
plot_label <- gsub('.F','', gsub('DPDP.','',rownames(samp_wk24)))
text(get_self(samp_wk24), get_best(samp_wk24) + 0.005, labels = plot_label)


lowcounts_samples <- c('DPDP.DO2.439.F','DPDP.DO2.473.F','DPDP.DO2.496.F','DPDP.DO2.572.F','DPDP.DO2.702.F','DPDP.DO2.719.F','DPDP.DO2.734.F','DPDP.DO2.727.F','DPDP.DO2.717.F')
points(get_self(samp_wk24)[lowcounts_samples], get_best(samp_wk24)[lowcounts_samples], pch = 21, bg = 'green')
mixture_samples   <- c('DPDP.DO2.549.F','DPDP.DO2.846.F')
points(get_self(samp_wk24)[mixture_samples], get_best(samp_wk24)[mixture_samples], pch = 21, bg = 'red')




