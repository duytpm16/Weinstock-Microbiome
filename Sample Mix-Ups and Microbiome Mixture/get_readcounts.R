options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)





### Comman line arguments / variables to change
args            <- commandArgs(trailingOnly = TRUE)
start           <- as.numeric(args[1])
end             <- as.numeric(args[2])
snp_list        <- readRDS('cc_variants_snp_list.rds')                               # Read in cc_variant list as generated from get_cc_variants.R
fastq_directory <- '/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/'  # Directory where fastq files are stored
chromosomes     <- c("1","2","3","4","5","6","7","8","9","10","11",                  # Vector of chromosomes
                     "12","13","14","15","16","17","18","19","X") 










### Change to directory where fastq files are stored
setwd(fastq_directory)
















### Loop through each MGS_DO_* directory
for(i in start:end){
    
    # Change to MGS_DO_* directory if it exists
    sample <- paste0('MGS_DO_',i)
    if(sample %in% list.dirs(fastq_directory, recursive = FALSE, full.names = FALSE)){
       samples_directory <- paste0(fastq_directory,sample)
       print(samples_directory)
       setwd(samples_directory)
    
    
    
    
    
      
      
       ### Loop through each week directory within the MGS_DO_* directory 
       for(week in list.dirs(samples_directory, recursive = FALSE, full.names = FALSE)){
           week_directory <- paste0(samples_directory,'/',week)
           print(week_directory)
           setwd(week_directory)
      
      
           pileup <- readRDS(grep('*pileup.rds', dir(), value = TRUE))
        
      
        
        
        
        
        
           ### Loop through each chromosome
           read_counts_list <- list()
           for(chr in chromosomes){
          
               # Get pileup and snpinfo of one chromosome
               chr_pileup     <- pileup[[chr]]
               chr_snp_info   <- snp_list[[chr]]
               chr_snp_info   <- chr_snp_info[ chr_snp_info$pos %in% chr_pileup$pos,]
            
            
            

              

    	       ### Create dataframe for major/minor allele
   	       major_read_counts_df <- chr_snp_info %>% 
                                                    select(pos, alleles) %>%
                                                    separate(col = alleles, into = c("nucleotide", "Minor"), sep = "\\|") %>%
                                                    select(-Minor)
               minor_read_counts_df <- chr_snp_info %>% 
                                                    select(pos, alleles) %>%
                                                    separate(col = alleles, into = c("Major", "nucleotide"), sep = "\\|") %>%
                                                    select(-Major)	          
            
               



               




                # Count number of reads that overlap snp
                snp_counts <- chr_pileup %>% 
                                         group_by(pos, nucleotide) %>%
                                         tally(count) %>%
                                         mutate(nucleotide = as.character(nucleotide))
  

    
                major_read_counts_df <- merge(major_read_counts_df, snp_counts, by = c('pos','nucleotide'), all.x = TRUE) %>%
                                                 dplyr::rename(Major_Allele = nucleotide, Major_Count = n)
                major_read_counts_df[is.na(major_read_counts_df)] <- 0
    

    

                minor_read_counts_df <- merge(minor_read_counts_df, snp_counts, by = c('pos','nucleotide'), all.x = TRUE) %>% 
                                                 dplyr::rename(Minor_Allele = nucleotide, Minor_Count = n)
                minor_read_counts_df[is.na(minor_read_counts_df)] <- 0
                  







                # Save read_counts_df into a list 
   	        read_counts_df <- merge(major_read_counts_df, minor_read_counts_df, by = 'pos')                
                read_counts_list[[chr]] <- read_counts_df

                
            } # For chr
           
                     
  
  
         saveRDS(read_counts_list,  file = paste0(week_directory, '/', sample, '_', week, '_readcounts.rds'))     
      
      } # For week
       
   } # If
    
} # For sample


