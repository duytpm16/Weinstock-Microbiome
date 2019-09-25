options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)





### Comman line arguments / variables to change
args            <- commandArgs(trailingOnly = TRUE)
start           <- as.numeric(args[1])
end             <- as.numeric(args[2])
snp_dir         <- '/projects/churchill-lab/data/Weinstock/Pomp_Benson/fastq/'
fastq_directory <- '/projects/churchill-lab/data/Weinstock/Pomp_Benson/fastq/'
chromosomes     <- as.character(1:19) 










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
           week_directory <- paste0(samples_directory,'/',week,'/bowtie1_run_v2')
           print(week_directory)
           setwd(week_directory)
      
      
      
      
      
           ### Loop through each chromosome
           for(chr in chromosomes){
        
               # Get pileup and snpinfo of one chromosome
               load(paste0(snp_dir, 'imp_snp_',chr,'.Rdata'))
               
               chr_pileup  <- readRDS(paste0(week_directory, '/', sample, '_', week, '_pileup_chr_',chr,'.rds'))
               snpinfo$pos <- round(snpinfo$pos * 1000000)
               snpinfo     <- snpinfo[snpinfo$pos %in% chr_pileup$pos,]      
    
        
        
        
        
        
        
               # Count number of reads that overlap snp
               snp_counts <- chr_pileup %>% 
                                group_by(pos, nucleotide) %>%
                                tally(count) %>%
                                mutate(nucleotide = as.character(nucleotide))
        
        
        
               allele_1 <- merge(snpinfo, snp_counts, by.x = c('pos','allele1'), by.y = c('pos', 'nucleotide'), all.x = TRUE) 
               allele_1 <- allele_1 %>% rename(n1 = n)
               allele_1$n1[is.na(allele_1$n1)] <- 0
               
               allele_2 <- merge(snpinfo, snp_counts, by.x = c('pos','allele2'), by.y = c('pos', 'nucleotide'), all.x = TRUE) 
               allele_2 <- allele_2 %>% rename(n2 = n)
               allele_2$n2[is.na(allele_2$n2)] <- 0      
        
  
        
        
               # Merge the two counts as one 
               counts <- merge(allele_1[,c('chr','pos','allele1','n1')], allele_2[,c('chr','pos','allele2','n2')], by = c('chr','pos'))                
               counts <- counts %>% arrange(pos)
               stopifnot(all(counts$pos %in% snpinfo$pos))
        
        
               # Save
               saveRDS(counts,  file = paste0(week_directory, '/', sample, '_', week, '_readcounts_chr_',chr,'.rds'))
        
        
           } # For chr
      
        } # For week
    
    } # If

} # For sample
