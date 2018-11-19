#####################################################################################################################
#   
#    This script creates directories to help organize the fastq files given by Jethro via Globus.
#    Note*: According to Jethro, 
#               1.) fastq.1 and fastq.2 are paired-end reads
#               2.) fastq.3 are singletons whose mate pair was lost during trimming                    
#   
#
#    Input:
#       1.) Path to where directories of each sample was created.
#
#
#    Output:
#       1.) Sub-directories within each of the sample's directory. 
#           Just moving files around to help organize the whole host_fastq directory.
#
#
#
#    Author: Duy Pham
#    E-mail: duy.pham@jax.org
#    Date:   November 14, 2018
#
#####################################################################################################################
library(filesstrings)



### Variables to change
orig_dir <- '/projects/churchill-lab/data/Weinstock/Pomp_Benson/host_fastq/'
setwd(orig_dir)






### Within each of the sample's directory, create a directory for week 6, 17, and 24, and then within each of the week's
#       directory, create a Pair_End and Singleton directory. Move files to their appropriate directory.
for(i in dir()){


    # Set directory to one of the MGS_DO_* directory
    mwgs_directory <- paste0(orig_dir,'/',i,'/')
    setwd(mwgs_directory)





    # Make a separate directory for each of the three weeks
    dir.create(paste0(mwgs_directory,'Week_6/'))
    dir.create(paste0(mwgs_directory,'Week_17/'))
    dir.create(paste0(mwgs_directory,'Week_24/'))


   
    
    
    # Moving files to their corresponding week directory
    week_6_mwgs <- dir()[grep('*-w6-*',dir())]
    week_17_mwgs <- dir()[grep('*-w17-*',dir())]    
    week_24_mwgs <- dir()[grep('*-w24-*',dir())]
    
    sapply(week_6_mwgs, FUN = function(x) file.move(x, paste0(mwgs_directory,'Week_6/')))
    sapply(week_17_mwgs, FUN = function(x) file.move(x, paste0(mwgs_directory,'Week_17/'))) 
    sapply(week_24_mwgs, FUN = function(x) file.move(x, paste0(mwgs_directory,'Week_24/')))






    ### For each of the week directory, create a Pair_End and Singleton directory and move files to their corresponding directory

    #   According to Jethro:
    #      fastq.1 and fastq.2 files are pair-end reads. Corresponding reads appear at the same position in each file
    #      fastq.3 are singleton reads whose mate pair was lost during trimming


    #   Do above for week 6 fastq files
    setwd(paste0(mwgs_directory,'/Week_6/'))                                                # Change directory to the week 6 directory
 
    dir.create(paste0(getwd(),'/Pair_End/'))                                                # Create Pair_End directory
    dir.create(paste0(getwd(),'/Singletons/'))                                              # Create singletons directory
 
    pair_end_fastq <- dir()[grep('*fastq.[12]', dir())]                                     # Find files that end in fastq.1 and fastq.2 for pair-end reads
    singleton_fastq <- dir()[grep('*fastq.3', dir())]                                       # Find files that end in fastq.3 for singletons

    sapply(pair_end_fastq, FUN = function(x) file.move(x, paste0(getwd(),'/Pair_End/')))    # Move pair-end reads to Pair_End directory
    file.move(singleton_fastq, paste0(getwd(),'/Singletons/'))                              # Move singleton reads to Singletons directory
 
    setwd(mwgs_directory)                                                                   # Move back to the MGS_DO_* directory
     
       	      										    # Above lines are repeated for week 17 and 24
    
 

    #   Do above for week 17 fastq files
    setwd(paste0(mwgs_directory,'/Week_17/'))
    
    dir.create(paste0(getwd(),'/Pair_End/'))
    dir.create(paste0(getwd(),'/Singletons/')) 
   
    pair_end_fastq <- dir()[grep('*fastq.[12]', dir())]
    singleton_fastq <- dir()[grep('*fastq.3', dir())]
    
    sapply(pair_end_fastq, FUN = function(x) file.move(x, paste0(getwd(),'/Pair_End/')))
    file.move(singleton_fastq, paste0(getwd(),'/Singletons/'))
    
    setwd(mwgs_directory)

    



    #  Do above for week 24 fastq files
    setwd(paste0(mwgs_directory,'/Week_24/'))
    
    dir.create(paste0(getwd(),'/Pair_End/'))
    dir.create(paste0(getwd(),'/Singletons/')) 
   
    pair_end_fastq <- dir()[grep('*fastq.[12]', dir())]
    singleton_fastq <- dir()[grep('*fastq.3', dir())] 

    sapply(pair_end_fastq, FUN = function(x) file.move(x, paste0(getwd(),'/Pair_End/')))
    file.move(singleton_fastq, paste0(getwd(),'/Singletons/'))
}
        

