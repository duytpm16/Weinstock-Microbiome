### Options and libraries
options(stringsAsFactors = FALSE)
library(RSQLite)
library(dplyr)
library(qtl2)







### Load 36 state genoprobs and sqlite cc variant file. Specify chromosomes
load('weinstock_pomp_36_state_genoprob.Rdata')
db  <- dbConnect(SQLite(), 'cc_variants.sqlite')
chr <- c('1',2:19)









### Imputed SNPs begin
for(i in chr){
  
    # Get SNPs for chromosome
    snpinfo <- dbGetQuery(db, paste0("SELECT * FROM variants WHERE chr=='", i, "'"))
    snpinfo$pos <- snpinfo$pos / 1e6
    
    
    # Remove complex SNPs (Ex. 'A|G/C')
    complex <- grepl("/", snpinfo$alleles, fixed = TRUE)
    snpinfo <- snpinfo[!complex,]
    
    
  
  
  
    
    
    # Imput snps for chromosome. Need to break into chunks to run if there are over 1m SNPs
    if(nrow(snpinfo) > 1e6){
       snps <- split(snpinfo, (as.numeric(rownames(snpinfo)) - 1) %/% 1e6)
       
       for(j in 1:length(snps)){
           index_snp <- snps[[j]] %>% index_snps(map = map, snpinfo = .) %>% mutate(index = 1:nrow(.))
           snp_pr    <- genoprob_to_snpprob(genoprobs = genoprob, snpinfo = index_snp)
           
           if(j == 1){
              imp_snps <- maxmarg(probs = snp_pr, cores = 32)[[1]]
           }else{
              imp_snps <- cbind(imp_snps, maxmarg(probs = snp_pr, cores = 32)[[1]])
           }
       }
       
      
     
    }else{
       index_snp <- snpinfo %>% index_snps(map = map, snpinfo = .) %>% mutate(index = 1:nrow(.)) 
       snp_pr    <- genoprob_to_snpprob(genoprobs = genoprob, snpinfo = index_snp)
       imp_snps  <- maxmarg(snp_pr, cores = 32)[[1]]
    }
  
    print(dim(imp_snps))

    
    
  
  
  
    
    
    # Pull the alleles from the SNP table
    alleles <- strsplit(snpinfo$alleles, "\\|")
    snpinfo$allele1 <- sapply(alleles, "[", 1)
    snpinfo$allele2 <- sapply(alleles, "[", 2)
    
    
    # Save as .Rdata file for chromosome
    save(imp_snps, snpinfo, file=paste0("/projects/churchill-lab/data/Weinstock/Pomp_Benson/fastq/imp_snp_", i, ".Rdata"))
    
    print(i)
}
