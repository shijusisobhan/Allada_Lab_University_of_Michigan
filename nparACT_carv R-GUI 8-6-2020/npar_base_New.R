
#***************************************************************************************
# This is the modified version of npar_base. This version retun not only nparact result 
#but also data for ploting

#Eg-:  npar_base_New("sleepstudy", 1/60, cutoff = 1)
# where sleepstudy is the actgraphy data with dat-time string and activit

#***************************************************************************************


npar_base_New<-function(name, SR, cutoff = 1){


  
data <- get(name)
if (is.data.frame(data)==F){
  data = as.data.frame(data)
}
if(ncol(data) == 2){
  data[,1] <- as.POSIXct(data[,1])
  data[,2] <- as.numeric(as.character(data[,2]))
  names(data)[1] <- "time"
  names(data)[2] <- "activity"
} 
if(ncol(data) == 3){
  names(data)[1] <- "date"
  names(data)[2] <- "time"
  names(data)[3] <- "activity"
  data$date <- NULL
  data$time <- as.POSIXct(data$time, format="%H:%M:%S")  
  data$activity <- as.numeric(as.character(data$activity))
}
if (any(is.na(data$activity)) == TRUE) stop("Please check your data! It must not contain NAs")

bin_hr <- 60  
a <- nrow(data) 
e <- SR*60 ## samples per minute
m <- bin_hr*SR*60  ## samples per hour
full_days <- floor(a/(e*bin_hr*24))

## --- Cut data to full days

  data <- data[1:(e*bin_hr*24*full_days),]

a <- nrow(data) 
b <- floor(a/(SR*60)) ## full minutes recorded


nparACT_auxfunctions1$nparACT_filt(data, a, cutoff)

if (SR != 1/60){
  data_min <- nparACT_auxfunctions1$nparACT_data_min(b, SR, data)
}  else {
  data_min <- data$activity
}


data_hrs <- nparACT_auxfunctions1$nparACT_data_hrs(data, a, m)
## -----------------------------------------------------------------------------


result_ISIV <- nparACT_ISIVfunctions$nparACT_ISIV(data_hrs, bin_hr)
IS <- result_ISIV[1]
IV <- result_ISIV[2]



minaverage <- nparACT_auxfunctions1$nparACT_minaverage(a, data_min)


## ---- L5, M10, RA calculation
result_RA <- nparACT_RAfunctions$nparACT_L5M10(data, minaverage, a, SR)
L5 <- result_RA[1]
L5_starttime <- result_RA[2]
M10 <- result_RA[3]
M10_starttime <- result_RA[4]
RA <- result_RA[5]

nparACT_result <- data.frame(IS, IV, RA, L5, L5_starttime, M10, M10_starttime)


nparACT_result_all <- list(data,minaverage,a,nparACT_result)

## --------------------------------

}



