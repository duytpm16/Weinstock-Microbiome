library(qtl2)
library(dplyr)



### Variables to change
args = commandArgs(trailing = TRUE)
query_variants <- create_variant_query_func("/projects/churchill-lab/resource/qtl2_snps_genes/cc_variants.sqlite", filter = "type=='snp'")
genoprobs      <- readRDS("pomp_genoprobs_qtl2.rds")
maps           <- readRDS("pomp_map_qtl2.rds")
chromsomes     <- c(args[1])






### SNP imputation begin
for(i in chromsomes){
    
    
    # Names to assign snp info and imputed snp to global environment
    snpinfo_name  <- paste0('snpinfo_chr', i)
    imp_snps_name <- paste0('imp_snps_chr', i)
  
  
  
    
    
    
    # Get snp information of one chromosome
    snpinfo_chr <- query_variants(chr = i, start = 1, end = 500)
    complex     <- grepl("/", snpinfo_chr$alleles, fixed=TRUE)
    snpinfo_chr <- snpinfo_chr[!complex,]
    assign(snpinfo_name, snpinfo_chr)
    

    
    
    
    
    
    # Indexing SNP information
    indexed_snpinfo <-  snpinfo_chr %>%                              
                                    rename(snp = snp_id)        %>%   # Renaming snp_id column to snp for the function ?index_snps
                                    index_snps(map = maps,            # Indexing snps (qtl2)
                                               snpinfo = .)      %>%    
                                    mutate(index = 1:nrow(.))         # Forcing next codes to work for all snps
    
             
  
    
    
    
    ### genoprob_to_snpprob doesn't seem to work when snpinfo matrix contains over 2m rows. So I am breaking it into two dataframes
    n <- nrow(indexed_snpinfo)
    if(n > 1600000){
       
      
       # Breaking snpinfo into two dataframes and then re-indexing
       indexed_snpinfo_1 <- snpinfo_chr[1:1600000,] %>%                              
                                                    index_snps(map = maps,      
                                                               snpinfo = .)     %>% 
                                                    mutate(index = 1:nrow(.))       
      
       indexed_snpinfo_2 <- snpinfo_chr[1600001:n,] %>%                              
                                                    index_snps(map = maps,          
                                                               snpinfo = .)    %>%    
                                                    mutate(index = 1:nrow(.))       
      
       
       
       
       # Compute SNP probabilities
       snp_pr1 <- genoprob_to_snpprob(genoprobs = genoprobs, snpinfo = indexed_snpinfo_1)
       snp_pr2 <- genoprob_to_snpprob(genoprobs = genoprobs, snpinfo = indexed_snpinfo_2)
       
       
       # Impute SNPs. Drop list and get the matrix '[[1]]' 
       imp_snps1 <- maxmarg(snp_pr1, cores = 0)[[1]]
       imp_snps2 <- maxmarg(snp_pr2, cores = 0)[[1]]
       
       
       # Making sure rownames match
       stopifnot(rownames(imp_snps1) == rownames(imp_snps2))
       
       # Cbind imputed snp matrices
       assign(imp_snps_name, cbind(imp_snps1, imp_snps2))
       
       
       
       
       
    }else{
      
       # Compute SNP probabilities
       snp_pr <- genoprob_to_snpprob(genoprobs = genoprobs[,i], snpinfo = indexed_snpinfo)
       
      
       # Impute SNPs. Drop list and get the matrix '[[1]]' 
       assign(imp_snps_name, maxmarg(snp_pr, cores = 0)[[1]])
       
    }

    
  
    
    

    
    
    
    
    
    
    # Save imputed snp matrix and snp information matrix as RData
    save(list = c(imp_snps_name, snpinfo_name), file = paste0('imputed_snps_chr_',i,'.RData'))
    rm(list = c(imp_snps_name, snpinfo_name))
}




