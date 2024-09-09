
#************************************************************

# This function merges all the files those have same subject ID
# Eg:- readActiwareCSVs_action(dir(PATH, full.names=T))


#************************************************************



readActiwareCSVs_action <- function(filenames,idcodes=basename(filenames)){
  
 
  
  if(length(filenames)!=length(idcodes)){
    stop("Need an idcode for each file; filenames and idcodes should have same length")
  }
  cat("Reading",length(filenames),"files...\n")
  
  
  out <- lapply(1:length(filenames),function(i){Merge.ActiwareCSV(filenames[i])})
  
  out <- do.call(rbind,out)
  return(out)
}







