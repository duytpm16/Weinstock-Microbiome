### Options and libraries
options(stringsAsFactors = FALSE)
library(RSQLite)
library(dplyr)
library(qtl2)




### Command line argument: specify chromosome and cores
args  <- commandArgs(trailing = TRUE)
chrom <- args[1] 
cores <- as.numeric(args[2])





### Load 36 state genoprobs and sqlite cc variant file. Specify chromosomes
load('weinstock_pomp_36_state_genoprob.Rdata')
db <- dbConnect(SQLite(), '/projects/churchill-lab/resource/qtl2_snps_genes/cc_variants.sqlite')









# Get SNPs for chromosome
snpinfo <- dbGetQuery(db, paste0("SELECT * FROM variants WHERE chr=='", chrom, "'"))
snpinfo$pos <- snpinfo$pos / 1e6


# Remove complex SNPs (Ex. 'A|G/C')
complex <- grepl("/", snpinfo$alleles, fixed = TRUE)
snpinfo <- snpinfo[!complex,]










### Imput snps for chromosome. 
#    Need to break into chunks to run if there are over 1m SNPs
if(nrow(snpinfo) > 1e6){
  
   # Split snpinfo into batches of 1m
   snps <- split(snpinfo, (as.numeric(rownames(snpinfo))) %/% 1e6)
  
   # Index and compute snp probabilities for each batch
   for(j in 1:length(snps)){
       index_snp <- snps[[j]] %>% index_snps(map = map, snpinfo = .) %>% mutate(index = 1:nrow(.))
       snp_pr    <- genoprob_to_snpprob(genoprobs = genoprob, snpinfo = index_snp)
      
       # Turn list object into matrix
       snps[[j]] <- as.data.frame(maxmarg(probs = snp_pr, cores = cores)[[1]])
   }
  
   # Combind all results as 1 matrix
   imp_snps <- bind_cols(snps)
   rownames(imp_snps) <- rownames(snps[[1]])
   imp_snps <- as.matrix(imp_snps)
  
  
  
}else{
   index_snp <- snpinfo %>% index_snps(map = map, snpinfo = .) %>% mutate(index = 1:nrow(.)) 
   snp_pr    <- genoprob_to_snpprob(genoprobs = genoprob, snpinfo = index_snp)
   imp_snps  <- maxmarg(snp_pr, cores = cores)[[1]]
}


print(dim(imp_snps))










# Pull the alleles from the SNP table
alleles <- strsplit(snpinfo$alleles, "\\|")
snpinfo$allele1 <- sapply(alleles, "[", 1)
snpinfo$allele2 <- sapply(alleles, "[", 2)
snpinfo <- snpinfo[snpinfo$snp_id %in% colnames(imp_snps),]



# Save as .Rdata file for chromosome
save(imp_snps, snpinfo, file=paste0("/projects/churchill-lab/data/Weinstock/Pomp_Benson/genotypes/imputed_snps/imp_snp_", chrom, ".Rdata"))
