library(lineup2)
library(broman)



raw_counts <- readRDS('raw_counts.rds')
alignment_counts_wk6 <- readRDS('alignment_counts_week_6.rds')
alignment_counts_wk17 <- readRDS('alignment_counts_week_17.rds')
alignment_counts_wk24 <- readRDS('alignment_counts_week_24.rds')
samp_wk6   <- readRDS("Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk6.rds")
samp_wk17  <- readRDS("Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk17.rds")
samp_wk24  <- readRDS("Weinstock Microbiome/Mixups/single_results_prop_mismatch_wk24.rds")
pair_result_6  <- readRDS('mixture_results_wk_6.rds')
pair_result_17 <- readRDS('mixture_results_wk_17.rds')
pair_result_24 <- readRDS('mixture_results_wk_24.rds')
# samp_p_w6 <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk6_all_v2.rds")
# samp_p_w17  <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk17_all_v2.rds")
# samp_p_w24  <- readRDS("~/Desktop/Weinstock Microbiome/Mixups/sample_results_wk24_all_v2.rds")







### Week 6
grayplot(x = get_self(samp_wk6) + .001, y = get_best(samp_wk6) + .001,
         xlab="Proportion Discordant with Self", ylab="Minimum Proprotion Discordant",
         xlim=c(0, 0.213), ylim=c(0, 0.213), 
         xaxs="i", yaxs="i",
         main = 'Proportion of Discordant Week 6')


lowcounts_samples_w6 <- c('DPDP.DO2.410.F','DPDP.DO2.413.F','DPDP.DO2.429.F','DPDP.DO2.558.F','DPDP.DO2.567.F','DPDP.DO2.568.F','DPDP.DO2.569.F','DPDP.DO2.570.F','DPDP.DO2.571.F','DPDP.DO2.572.F','DPDP.DO2.806.F')
points(get_self(samp_wk6)[lowcounts_samples_w6] + .001, get_best(samp_wk6)[lowcounts_samples_w6] + .001, pch = 21, bg = 'green')
text(get_self(samp_wk6)[lowcounts_samples_w6[-c(2,9)]], get_best(samp_wk6)[lowcounts_samples_w6[-c(2,9)]] + 0.005, labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',lowcounts_samples_w6[-c(2,9)]))))
text(get_self(samp_wk6)[lowcounts_samples_w6[2]] - .01, get_best(samp_wk6)[lowcounts_samples_w6[2]] + 0.002, labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',lowcounts_samples_w6[2]))))
text(get_self(samp_wk6)[lowcounts_samples_w6[9]] + .012, get_best(samp_wk6)[lowcounts_samples_w6[9]] + 0.002, labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',lowcounts_samples_w6[9]))))
mixture_samples_w6   <- c('DPDP.DO2.444.F','DPDP.DO2.492.F','DPDP.DO2.541.F','DPDP.DO2.507.F','DPDP.DO2.552.F')
points(get_self(samp_wk6)[mixture_samples_w6] + .001, get_best(samp_wk6)[mixture_samples_w6] + .001, pch = 21, bg = 'red')
text(get_self(samp_wk6)[mixture_samples_w6] - 0.007, get_best(samp_wk6)[mixture_samples_w6] + 0.001 , labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',mixture_samples_w6))))
legend('topright', legend=c("Low Read Counts", "Mixtures", 'Okay'), col=c("green", "red", 'lightblue'), pch = 16, cex=1.5)








par(mfrow = c(4,2))
sample <- 'DPDP.DO2.444.F'
grayplot(pair_result_6[[sample]][,1], pair_result_6[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_6[[sample]][,1])     
max.y = which.max(pair_result_6[[sample]][,4])
val.x = pair_result_6[[sample]][max.x,1]
val.y = pair_result_6[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.04, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_6[[sample]])[[1]][max.x])))

plot(samp_p_w6[sample,], main = 'DO-444', pch = 21, col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w6[sample,][order(samp_p_w6[sample,])][1:5]
min_index  <- match(min_values, samp_p_w6[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk6)[min_index])))





sample <- 'DPDP.DO2.492.F'
grayplot(pair_result_6[[sample]][,1], pair_result_6[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-492',
         bgcolor = 'white')
max.x = which.max(pair_result_6[[sample]][,1])     
max.y = which.max(pair_result_6[[sample]][,4])
val.x = pair_result_6[[sample]][max.x,1]
val.y = pair_result_6[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.03, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_6[[sample]])[[1]][max.x])))

plot(samp_p_w6[sample,], main = 'DO-492', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w6[sample,][order(samp_p_w6[sample,])][1:5]
min_index  <- match(min_values, samp_p_w6[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk6)[min_index])))




sample <- 'DPDP.DO2.507.F'
grayplot(pair_result_6[[sample]][,1], pair_result_6[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-507',
         bgcolor = 'white')
max.x = which.max(pair_result_6[[sample]][,1])     
max.y = which.max(pair_result_6[[sample]][,4])
val.x = pair_result_6[[sample]][max.x,1]
val.y = pair_result_6[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.04, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_6[[sample]])[[1]][max.x])))


plot(samp_p_w6[sample,], main = 'DO-507', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w6[sample,][order(samp_p_w6[sample,])][1:5]
min_index  <- match(min_values, samp_p_w6[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk6)[min_index])))



sample <- 'DPDP.DO2.552.F'
grayplot(pair_result_6[[sample]][,1], pair_result_6[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-552',
         bgcolor = 'white')
max.x = which.max(pair_result_6[[sample]][,1])     
max.y = which.max(pair_result_6[[sample]][,4])
val.x = pair_result_6[[sample]][max.x,1]
val.y = pair_result_6[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.03, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_6[[sample]])[[1]][max.x])))



plot(samp_p_w6[sample,], main = 'DO-552', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w6[sample,][order(samp_p_w6[sample,])][1:5]
min_index  <- match(min_values, samp_p_w6[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk6)[min_index])))
par(mfrow = c(1,1))








### Week 17
grayplot(get_self(samp_wk17), get_best(samp_wk17),
         xlab="Proportion Discordant with Self", ylab="Minimum Proprotion Discordant",
         xlim=c(0, 0.213), ylim=c(0, 0.213), xaxs="i", yaxs="i",
         main = 'Proportion Discordant Week 17')


 

lowcounts_samples_w17 <- c('DPDP.DO2.415.F','DPDP.DO2.425.F','DPDP.DO2.741.F','DPDP.DO2.745.F','DPDP.DO2.550.F','DPDP.DO2.738.F','DPDP.DO2.739.F','DPDP.DO2.740.F','DPDP.DO2.742.F','DPDP.DO2.743.F')
points(get_self(samp_wk17)[lowcounts_samples_w17], get_best(samp_wk17)[lowcounts_samples_w17], pch = 21, bg = 'green')
text(get_self(samp_wk17)[lowcounts_samples_w17], get_best(samp_wk17)[lowcounts_samples_w17] + .003, labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',lowcounts_samples_w17))))
mixture_samples_w17   <- c('DPDP.DO2.440.F','DPDP.DO2.484.F', 'DPDP.DO2.833.F','DPDP.DO2.840.F','DPDP.DO2.841.F','DPDP.DO2.835.F','DPDP.DO2.836.F','DPDP.DO2.850.F')
points(get_self(samp_wk17)[mixture_samples_w17], get_best(samp_wk17)[mixture_samples_w17], pch = 21, bg = 'red')
text(get_self(samp_wk17)[mixture_samples_w17[-c(4,6,8)]] +.008, get_best(samp_wk17)[mixture_samples_w17[-c(4,6,8)]], labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',mixture_samples_w17[-c(4,6,8)]))))
text(get_self(samp_wk17)[mixture_samples_w17[c(4,6,8)]] -.008, get_best(samp_wk17)[mixture_samples_w17[c(4,6,8)]], labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',mixture_samples_w17[c(4,6,8)]))))
legend('topright', legend=c("Low Read Counts", "Mixtures", 'Okay'), col=c("green", "red",'lightblue'),pch = 16, cex=1.5)




par(mfrow = c(4,2))

sample <- 'DPDP.DO2.440.F'
grayplot(pair_result_17[[sample]][,1], pair_result_17[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_17[[sample]][,1])     
max.y = which.max(pair_result_17[[sample]][,4])
val.x = pair_result_17[[sample]][max.x,1]
val.y = pair_result_17[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.035, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_17[[sample]])[[1]][max.x])))

plot(samp_p_w17[sample,], main = 'DO-440', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w17[sample,][order(samp_p_w17[sample,])][1:5]
min_index  <- match(min_values, samp_p_w17[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk17)[min_index])))


sample <- 'DPDP.DO2.484.F'
grayplot(pair_result_17[[sample]][,1], pair_result_17[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_17[[sample]][,1])     
max.y = which.max(pair_result_17[[sample]][,4])
val.x = pair_result_17[[sample]][max.x,1]
val.y = pair_result_17[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.02, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_17[[sample]])[[1]][max.x])))

plot(samp_p_w17[sample,], main = 'DO-484', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w17[sample,][order(samp_p_w17[sample,])][1:5]
min_index  <- match(min_values, samp_p_w17[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk17)[min_index])))


sample <- 'DPDP.DO2.833.F'
grayplot(pair_result_17[[sample]][,1], pair_result_17[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_17[[sample]][,1])     
max.y = which.max(pair_result_17[[sample]][,4])
val.x = pair_result_17[[sample]][max.x,1]
val.y = pair_result_17[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.007, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_17[[sample]])[[1]][max.x])))

plot(samp_p_w17[sample,], main = 'DO-833', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w17[sample,][order(samp_p_w17[sample,])][1:5]
min_index  <- match(min_values, samp_p_w17[sample,])
text(min_index-30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk17)[min_index])))



sample <- 'DPDP.DO2.835.F'
grayplot(pair_result_17[[sample]][,1], pair_result_17[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_17[[sample]][,1])     
max.y = which.max(pair_result_17[[sample]][,4])
val.x = pair_result_17[[sample]][max.x,1]
val.y = pair_result_17[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.025, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_17[[sample]])[[1]][max.x])))

plot(samp_p_w17[sample,], main = 'DO-835', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w17[sample,][order(samp_p_w17[sample,])][1:5]
min_index  <- match(min_values, samp_p_w17[sample,])
text(min_index-30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk17)[min_index])))


sample <- 'DPDP.DO2.836.F'
grayplot(pair_result_17[[sample]][,1], pair_result_17[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_17[[sample]][,1])     
max.y = which.max(pair_result_17[[sample]][,4])
val.x = pair_result_17[[sample]][max.x,1]
val.y = pair_result_17[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.025, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_17[[sample]])[[1]][max.x])))

plot(samp_p_w17[sample,], main = 'DO-836', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w17[sample,][order(samp_p_w17[sample,])][1:5]
min_index  <- match(min_values, samp_p_w17[sample,])
text(min_index-30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk17)[min_index])))



sample <- 'DPDP.DO2.840.F'
grayplot(pair_result_17[[sample]][,1], pair_result_17[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_17[[sample]][,1])     
max.y = which.max(pair_result_17[[sample]][,4])
val.x = pair_result_17[[sample]][max.x,1]
val.y = pair_result_17[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.008, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_17[[sample]])[[1]][max.x])))

plot(samp_p_w17[sample,], main = 'DO-840', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w17[sample,][order(samp_p_w17[sample,])][1:5]
min_index  <- match(min_values, samp_p_w17[sample,])
text(min_index-30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk17)[min_index])))


sample <- 'DPDP.DO2.841.F'
grayplot(pair_result_17[[sample]][,1], pair_result_17[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_17[[sample]][,1])     
max.y = which.max(pair_result_17[[sample]][,4])
val.x = pair_result_17[[sample]][max.x,1]
val.y = pair_result_17[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.002, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_17[[sample]])[[1]][max.x])))

plot(samp_p_w17[sample,], main = 'DO-841', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w17[sample,][order(samp_p_w17[sample,])][1:5]
min_index  <- match(min_values, samp_p_w17[sample,])
text(min_index-30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk17)[min_index])))





sample <- 'DPDP.DO2.850.F'
grayplot(pair_result_17[[sample]][,1], pair_result_17[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_17[[sample]][,1])     
max.y = which.max(pair_result_17[[sample]][,4])
val.x = pair_result_17[[sample]][max.x,1]
val.y = pair_result_17[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.015, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_17[[sample]])[[1]][max.x])))

plot(samp_p_w17[sample,], main = 'DO-850', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w17[sample,][order(samp_p_w17[sample,])][1:5]
min_index  <- match(min_values, samp_p_w17[sample,])
text(min_index-30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk17)[min_index])))
par(c(1,1))








### Week 24
grayplot(get_self(samp_wk24), get_best(samp_wk24),
         xlab="Proportion Discordant with Self", ylab="Minimum Proprotion Discordant",
         xlim=c(0, 0.213), ylim=c(0, 0.213), xaxs="i", yaxs="i",
         main = 'Proportion Discordant Week 24')



lowcounts_samples_w24 <- c('DPDP.DO2.573.F','DPDP.DO2.574.F','DPDP.DO2.575.F','DPDP.DO2.439.F','DPDP.DO2.473.F','DPDP.DO2.496.F','DPDP.DO2.572.F','DPDP.DO2.702.F','DPDP.DO2.719.F','DPDP.DO2.734.F','DPDP.DO2.727.F','DPDP.DO2.717.F','DPDP.DO2.711.F')
points(get_self(samp_wk24)[lowcounts_samples_w24], get_best(samp_wk24)[lowcounts_samples_w24], pch = 21, bg = 'green')
text(get_self(samp_wk24)[lowcounts_samples_w24[-9]] - .008,get_best(samp_wk24)[lowcounts_samples_w24[-9]], labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',lowcounts_samples_w24[-9]))))
text(get_self(samp_wk24)[lowcounts_samples_w24[9]] + .009,get_best(samp_wk24)[lowcounts_samples_w24[9]], labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',lowcounts_samples_w24[9]))))
mixture_samples_w24   <- c('DPDP.DO2.411.F','DPDP.DO2.416.F','DPDP.DO2.539.F','DPDP.DO2.549.F','DPDP.DO2.472.F','DPDP.DO2.846.F','DPDP.DO2.801.F','DPDP.DO2.732.F','DPDP.DO2.825.F','DPDP.DO2.504.F')
points(get_self(samp_wk24)[mixture_samples_w24], get_best(samp_wk24)[mixture_samples_w24], pch = 21, bg = 'red')
text(get_self(samp_wk24)[mixture_samples_w24[-c(1,3,2,5,9)]], get_best(samp_wk24)[mixture_samples_w24[-c(1,3,2,5,9)]] + .003, labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',mixture_samples_w24[-c(1,3,2,5,9)]))))
text(get_self(samp_wk24)[mixture_samples_w24[c(3,2,5)]]-.008, get_best(samp_wk24)[mixture_samples_w24[c(3,2,5)]], labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',mixture_samples_w24[c(3,2,5)]))))
text(get_self(samp_wk24)[mixture_samples_w24[c(1,9)]]+.008, get_best(samp_wk24)[mixture_samples_w24[c(1,9)]], labels = gsub('2[.]','-',gsub('.F','', gsub('DPDP.','',mixture_samples_w24[c(1,9)]))))
legend('topright', legend=c("Low Read Counts", "Mixtures", "Okay"), col=c("green", "red", "lightblue"), pch = 16, cex=1.5)




par(mfrow = c(4,2))
sample <- 'DPDP.DO2.411.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = gsub('.F','',gsub('2[.]','-',gsub('DPDP.','',sample))),
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.015, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))

plot(samp_p_w24[sample,], main = 'DO-411', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))





sample <- 'DPDP.DO2.416.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-416',
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.006, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))

plot(samp_p_w24[sample,], main = 'DO-416', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))






sample <- 'DPDP.DO2.472.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-472',
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.005, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))

plot(samp_p_w24[sample,], main = 'DO-472', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))





sample <- 'DPDP.DO2.539.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-539',
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.01, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))

plot(samp_p_w24[sample,], main = 'DO-539', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))




sample <- 'DPDP.DO2.549.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-549',
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.03, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))

plot(samp_p_w24[sample,], main = 'DO-549', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))




sample <- 'DPDP.DO2.732.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-732',
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.02, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))

plot(samp_p_w24[sample,], main = 'DO-732', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))




sample <- 'DPDP.DO2.801.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-801',
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.015, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))

plot(samp_p_w24[sample,], main = 'DO-801', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index+30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))






sample <- 'DPDP.DO2.846.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2-846',
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.04, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))


plot(samp_p_w24[sample,], main = 'DO-846', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index-30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))


sample <- 'DPDP.DO2.825.F'
grayplot(pair_result_24[[sample]][,1], pair_result_24[[sample]][,4],
         xlab=expression(hat(p)), 
         ylab="LRT for p=0", 
         main = 'DO2.504',
         bgcolor = 'white')
max.x = which.max(pair_result_24[[sample]][,1])     
max.y = which.max(pair_result_24[[sample]][,4])
val.x = pair_result_24[[sample]][max.x,1]
val.y = pair_result_24[[sample]][max.y,4]
points(val.x, val.y)
text(val.x-.04, val.y , labels = gsub('.F','',gsub('DPDP.','', dimnames(pair_result_24[[sample]])[[1]][max.x])))


plot(samp_p_w24[sample,], main = 'DO-846', pch = 21,  col = 'black', bg = 'lightblue', xlab = 'Genomic DNA Sample', ylab = 'Proportion Discordant')
min_values <- samp_p_w24[sample,][order(samp_p_w24[sample,])][1:5]
min_index  <- match(min_values, samp_p_w24[sample,])
text(min_index-30, min_values, labels = gsub('.[mMfF]','',gsub('DPDP.','',colnames(samp_wk24)[min_index])))

par(mfrow = c(1,1))












