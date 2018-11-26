library(Rsamtools)
library(qtl2)





### Read in cc variants list as generated from get_cc_variants.R
snp_list <- readRDS('cc_variants_snp_list.rds')






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
