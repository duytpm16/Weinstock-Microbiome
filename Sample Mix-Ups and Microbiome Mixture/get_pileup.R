library(Rsamtools)
library(qtl2)




### Getting the CC variants of each chromosome and storing to a list
query_variants <- create_variant_query_func("cc_variants.sqlite")

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







### Change to directory where fastq files are stored
fastq_directory <- '/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/'
setwd(fastq_directory)








### Loop through each MGS_DO_* directory
for(sample in dir()){
  
    samples_directory <- paste0(fastq_directory,sample)
    print(samples_directory)
    setwd(samples_directory)
  
     
     
     
     
     
    ### Loop through each week directory within the MGS_DO_* directory 
    for(week in grep('Week_[0-9]$', dir(), value = TRUE)){
      
        week_directory <- paste0(samples_directory, '/', week)
        print(week_directory)
        setwd(week_directory)
        
        
        # Get bam and indexed bam file
        bf <- BamFile(file  = grep('*.bam$', dir(), value = TRUE),
                      index = grep('*.bai$', dir(), value = TRUE))
        
        
        
        
      
        
  
        
        ### Loop through each chromosome to get pileup
        pileup_list <- list()
        for(chromosome in c('1')){
            chr <- paste0('chr',chromosome)                                                         # Need to paste 'chr' to chromosome since UCSC named their contigs with 'chr*'
            
            
            chr_length <- as.data.frame(seqinfo(bf))[chr, "seqlengths"]
            chr_param  <- ScanBamParam(which=setNames(IRangesList(IRanges(0L, chr_length)), chr))
            
            
            chr_pileup <- pileup(bf, scanBamParam=param)
            chr_pileup <- chr_pileup[chr_pileup$pos %in% snp_list[[i]]$pos, ]
            
            
            pileup_list[[chromosome]] <- chr_pileup
            
        } # For chromosome
        
        
        
        saveRDS(pileup_list, file = paste0(week_directory, '/', sample, '_', week, '_pileup.rds'))
    } # For week
   
  
  
} # For sample
