### Command line arguements / Variables to change
#   1.) results_dir - path to were the pair_sample_summary.R or single_sample_summary.R results are stored
#   2.) type - paired or sample
#   3.) week - 6, 17, or 24
args = commandArgs(trailingOnly = TRUE)
results_dir <- args[1]
type        <- args[2]
week        <- args[3]






### Change to directory where results are stored
setwd(results_dir)



combined_results <- NULL



for(chr in 1:19){
    
    # Get file name
    results_file <- paste0(type, '_results_week_', week, '_chr_', chr, ".rds")
    
    print(results_file)
    
    # Check if file exists
    if(file.exists(results_file)){
      
       # Read in file
       result <- readRDS(results_file) 
       
       
       
       
       
       # If no results have been stored to combined_results, store result, else combine result with combined_results 
       if(is.null(combined_results)){
      
          combined_results <- result
          
          
          
      
       }else{
          
          stopifnot(length(result) == length(combined_results))
         
          result <- result[names(combined_results)]
          
          for(i in seq_along(combined_results)) {
            
              stopifnot(dim(result[[i]]) == dim(combined_results[[i]]),
                        all(rownames(result[[i]]) == rownames(combined_results[[i]])))
            
            
              combined_results[[i]] <- combined_results[[i]] + result[[i]]
          }
       }
    }
  

}


saveRDS(combined_results, file=paste0(type,"_results_all.rds"))

