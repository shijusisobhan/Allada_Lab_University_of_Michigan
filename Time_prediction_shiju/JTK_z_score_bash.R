
rm(list=ls())

source("/home/shijusis/JTK_analysis/JTK_CYCLEv3.1.R")

data_all<-read.csv('/home/shijusis/JTK_analysis/Moller_archer_Broun_data.csv')


data_all<-data_all[-c(1,2),]
data_all<-na.omit(data_all[!(is.na(data_all) | data_all==""), ]) # remove blank and NAN

limt<-(ncol(data_all)-1)/2
data_all <- data_all[rowSums(data_all != 0) > limt,] # Eliminate row with 50% entries are 0

# data_all<-data_all[!duplicated(data_all$samp), ] # Remove duplicate names
rownames(data_all)<- NULL # Reset the row number

data=data_all[-1]
Gene<-data_all[1]



# Function to calculate z-scores for a row
z_score <- function(x) {
  (x - mean(x)) / sd(x)
}


# Apply the function to each row
data_z_scores <- t(apply(data, 1, z_score))

## *********************************************************************************

# Start JTK

# source("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/JTKversion3/JTK_CYCLEv3.1.R")

Time_rep<-c(29, 25, 32, 24, 23, 31, 25, 27, 31, 28, 29, 36, 28, 28,
            39, 33, 32, 46, 37, 34, 44, 37, 30, 38)

jtkdist(24, Time_rep)       # 12 total time points, replica given in Time_rep
periods <- 23:24       # looking for rhythms between 21-23 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,1)  # 2 is the number of hours between time points



flush.console()


st <- system.time({
  res <- apply(data_z_scores,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- res[order(res$ADJ.P,-res$AMP),]
})


res<-cbind(Gene, res)

write.csv(res, '/home/shijusis/JTK_analysis/JTK_z_score_all.csv')