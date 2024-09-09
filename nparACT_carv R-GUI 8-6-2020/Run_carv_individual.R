 

Run_carv_individual<-function(PATH, File_names, ind_st){
  
  # This is the function which run carv, npar analysis for each file (Individual file run)
  
source("carv.R")

ind_st.time<<-ind_st

my.CSV.files <- PATH



## Carv Run
actiware.data <<- readActiwareCSVs(my.CSV.files)
carv_individual<<-NLSredux(actiware.data)
carv_individual$idcode<-File_names

sub_split<-split(actiware.data, actiware.data$idcode)



Ind.samp.time=c()
for (isub in 1:length(sub_split)){
  
  Ind.samp.time[isub]<-(sub_split[[isub]]$time[2]-sub_split[[isub]]$time[1])*60
  
}

carv_individual$T.sampling<-Ind.samp.time

carv_individual<-carv_individual[,c("idcode","T.sampling", "missing.data", "actamp","actbeta","actphi","actmin","actmesor","actupmesor","actdownmesor","actalph","actwidthratio","RsqNLS","Fact","Fnlrgact")]



source("nparactiware.R")

outDir1= "npar_processed_files"

unlink("npar_processed_files", recursive = TRUE)  # delete previously created folder

SRR<<-round(1/(as.numeric(Ind.samp.time[1])*60), digits=5)

npar_result_individual<-mod_nparACT_base_loop(my.CSV.files, outDir1, SR=SRR, cutoff = 1, plot = T, fulldays = T)

unlink("npar_processed_files", recursive = TRUE)  # delete previously created folder

return(list(carv_individual, npar_result_individual))


}