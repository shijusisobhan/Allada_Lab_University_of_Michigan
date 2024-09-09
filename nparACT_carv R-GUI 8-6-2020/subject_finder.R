
#***********************************************************************
# This is the function to identfy the subjects of the files that uploaded
# Eg: subject_finder(dir(PATH, full.names=T))
#**********************************************************************

subject_finder <- function(filenames){
  
  
  
  Id_data=c()
  
  idrow1 <- lapply(1:length(filenames),function(i){
    
    tmp <- readLines(filenames[i])
    
    idRow <- grep("Identity:",tmp)
    
    Id_data1<- read.csv(filenames[i], skip=  idRow-1, nrows = 1, header = F)
    
    Id_data[i] =Id_data1[2]
    
  })
  

  idrow2=unique(idrow1)
  
  return(idrow2)
  
}