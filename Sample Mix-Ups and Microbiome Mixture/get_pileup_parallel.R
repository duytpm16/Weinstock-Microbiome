############################################################################################################################################################
#   
#    This script allows generating pile-up files in parallel for each sample directory
#
#       Note*: Directory are named like MGS_DO_*, where * is a number between 401-850.
#              Within each of the directory, there are 3 sub-directories: Week_6, Week_17, and Week_24.
#              Within each of the sub-directories, there is a .bam and .bai file.
#
#
#
#
#    Input:
#        1.) Start range of *. 
#        2.) End range of *.
#
#
#    Output:
#        1.) Pile-up files for each of the 3 weeks in the MGS_DO_* directory
#
#
#    Author: Duy Pham
#    E-mail: duy.pham@jax.org
#    Date  : November 25, 2018
#
############################################################################################################################################################
library(Rsamtools)
library(qtl2)


### Comman line arguments/variables to change
args = commandArgs(trailingOnly = TRUE)
start           <- as.numeric(args[1])
end             <- as.numeric(args[2])
snp_dir         <- '/home/phamd/'
fastq_directory <- '/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/'      # Directory where fastq files are stored






setwd(fastq_directory)









### Loop through each MGS_DO_* directory
for(i in start:end){

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
        
           # Get bam and indexed bam file
           bf <- BamFile(file  = grep('*.bam$', dir(), value = TRUE),
                         index = grep('*.bai$', dir(), value = TRUE))
        
        
        
        
        
  
        
           ### Loop through each chromosome to get pileup
           for(chromosome in c(1:19,'X')){
               load(paste0(snp_dir,'imputed_snps_chr_',chromsome,'.RData'))
               snpinfo     <- get(paste0('snpinfo_chr',chromsome))
               snpinfo$pos <- round(snpinfo$pos * 1000000)
               
               
               chr <- paste0('chr',chromosome)                                                         # Need to paste 'chr' to chromosome since UCSC named their contigs with 'chr*'
               
            
               chr_length <- as.data.frame(seqinfo(bf))[chr, "seqlengths"]
               chr_param  <- ScanBamParam(which=setNames(IRangesList(IRanges(0L, chr_length)), chr))
            
            
               chr_pileup <- pileup(bf, scanBamParam=chr_param)
               chr_pileup <- chr_pileup[chr_pileup$pos %in% snp_list[[chromosome]]$pos, ]
            
            
               saveRDS(pileup_list, file = paste0(week_directory, '/', sample, '_', week, '_pileup_chr_',chromsome,'.rds'))   
            
           } # For chromosome
        
       } # For week
   
    } # If sample exists

} # For i
