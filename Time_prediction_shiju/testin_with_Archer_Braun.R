

# before running this program run program_shiju.R first

#Data1<-read.csv('Braun_data_and_DLMO_all.csv')
Data1<-read.csv('Archer_data_and_DLMO_all.csv')
DLMO25<-t(Data1[1,-1])

cycling.expr<-Data1[which(Data1$ext_gene %in% cycling_genes$ext_genes ),-1]




## Z-Score
sample.mean <- apply(cycling.expr,2,mean)
sample.std <- apply(cycling.expr,2,sd)
z.expr <- apply(cycling.expr, 1, function(x) {
  tmp <- (x-sample.mean)/sample.std
})



## Compute ratio of genes or difference of log2 value
ratio.expr <- apply(cycling.expr,2,function(xx){
  tmp <- outer(xx,xx,'-')
  tmp[upper.tri(tmp)]
})



ratio.label <- apply(cycling.expr,2,function(xx){
  tmp <- outer(rownames(cycling.expr),rownames(cycling.expr),"paste")
  tmp[upper.tri(tmp)]
})

ratio.label <- apply(cycling.expr,2,function(xx){
  tmp <- outer(rownames(cycling.expr),rownames(cycling.expr),"paste")
  tmp[upper.tri(tmp)]
})

rownames(ratio.expr) <- ratio.label[,1]


ratio.expr<-t(ratio.expr) # genes on the column and ratios in the rows

#Conver time into angles --->   a=2*pi*time/24

DLMO25_angle<-(DLMO25%%24)*2*pi/24

# Convert angle into cartesian coordinates   
y1<-sin(DLMO25_angle) # sin(a)
y2<-cos(DLMO25_angle) # cos(a)
y=data.frame(y1,y2) # y= [sin(a) cos(a)]

Pred_net<-(predict(Ridge_regg, s=Ridge_regg$lambda.1se, newx = ratio.expr))
Pred_angle<-Pred_net[,,1] # Extract the predicted angle
Pred_Time<-(atan2(Pred_angle[,1],Pred_angle[,2])%%(2*pi))*(24/(2*pi)) # convert back angle to time
y_test_Time<-(atan2(y[,1],y[,2])%%(2*pi))*(24/(2*pi))


## ****plot original time vs predicted time ************

plot(y_test_Time, Pred_Time, main='Archer et al.',
     xlab='True', ylab='Predicted')


hrerr <- abs(Pred_Time-y_test_Time)
hrerr <- hrerr[!is.na(hrerr)]
hrerr <- pmin(hrerr,24-hrerr)
mederr=median(abs(hrerr))
hrsoff <- seq(0,12,length=49)



fracacc <- sapply(hrsoff,function(hrtol){
  100*sum(abs(hrerr)>hrtol)/length(hrerr)
})


norm.fracacc <- (100-fracacc)/100
norm.hrsoff <- hrsoff/12
auc <- round(sum(norm.fracacc[-1]*diff(norm.hrsoff)), digit=2)



plot(hrsoff,100-fracacc,type='l', xlim=c(0,12), xlab = "correct to within (hrs)",
     ylab= '% correct', main=  paste0("Archer et al. (", 'auc:', auc, ")"))
abline(a=0,b=100/12,col="grey")
