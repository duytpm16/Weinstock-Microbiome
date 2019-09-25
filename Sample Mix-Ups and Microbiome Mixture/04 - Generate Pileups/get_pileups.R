############################################################################################################################################################
#   
#    This script generate pile-up files in parallel for each DO sample
#
#       Note*: Directory are named like MGS_DO_*, where * is a number between 401-850.
#              Within each of the directory, there are 3 sub-directories: Week_6, Week_17, and Week_24.
#              Within each of the sub-directories, there is a bowtie1_run sub-directory that contains the sorted.bam and sorted.bam.bai files.
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
#        1.) Pile-up files for each chromsome at each week
#
#
#    Author: Duy Pham
#    E-mail: duy.pham@jax.org
#    Date  : November 25, 2018
#
############################################################################################################################################################
library(Rsamtools)



### Comman line arguments/variables to change
args = commandArgs(trailingOnly = TRUE)
start           <- as.numeric(args[1])
end             <- as.numeric(args[2])
snp_dir         <- '/projects/churchill-lab/data/Weinstock/Pomp_Benson/fastq/'
fastq_directory <- '/projects/churchill-lab/data/Weinstock/Pomp_Benson/fastq/'






setwd(fastq_directory)









### Loop through each MGS_DO_* directory
for(i in start:end){
  
    sample <- paste0('MGS_DO_',i)
    if(sample %in% list.dirs(fastq_directory, recursive = FALSE, full.names = FALSE)){
       samples_directory <- paste0(fastq_directory,sample)
       print(samples_directory)
       setwd(samples_directory)
    
    
    
    
    
       # Loop through each week directory within the MGS_DO_* directory 
       for(week in list.dirs(samples_directory, recursive = FALSE, full.names = FALSE)){
           week_directory <- paste0(samples_directory,'/',week,'/bowtie1_run_v2')
           print(week_directory)
           setwd(week_directory)
        
           # Get sorted bam and indexed file
           bf <- BamFile(file  = grep('*sorted.bam$', dir(), value = TRUE),
                         index = grep('*sorted.bam.bai$', dir(), value = TRUE))
        
        
        
        
        
        
        
           # Loop through each chromosome to get pileup
           for(chromosome in c(1:19)){
             
               # Load impute snps data
               load(paste0(snp_dir,'imp_snp_',chromosome,'.Rdata'))                
               snpinfo$pos <- round(snpinfo$pos * 1000000)
               
               
          
               # Find pileup for chromosome
	       chr <- paste0('chr',chromosome)
               chr_length <- as.data.frame(seqinfo(bf))[chr, "seqlengths"]
               chr_param  <- ScanBamParam(which=setNames(IRangesList(IRanges(0L, chr_length)), chr))
          
               chr_pileup <- pileup(bf, scanBamParam=chr_param)
               chr_pileup <- chr_pileup[chr_pileup$pos %in% snpinfo$pos, ]
          
          
               
               # Save
               saveRDS(chr_pileup, file = paste0(week_directory, '/', sample, '_', week, '_pileup_chr_',chromosome,'.rds'))
          
	 }# For chromosome
      } # For week
   } # If sample exists
} # For i
