library(qtl2)



query_variants <- create_variant_query_func("/projects/churchill-lab/resource/qtl2_snps_genes/cc_variants.sqlite", filter = "type=='snp'")





snp_list <- list()
for(i in c(1:19,'X')){


    # Get CC variants at chromosome i
    snp_info     <- query_variants(chr = i, start = 1, end = 1000)
    snp_info$pos <- snp_info$pos * 1000000                           # Need to multiple pos by 1e6 to match with pileup position later in the script



    # Remove complex variants
    complex  <- grepl("/", snp_info$alleles, fixed=TRUE)
    snp_info <- snp_info[!complex, ]



    snp_list[[i]] <- snp_info


}





saveRDS(snp_list, 'cc_variants_snp_list.rds')