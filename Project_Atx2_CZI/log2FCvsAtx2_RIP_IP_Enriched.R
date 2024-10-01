rm(list=ls())

DEG_Atx2<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Atx2/Writeup/Table-2.csv')

Intersect_gens<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Atx2/Writeup/Atx2_Wake_sleep/Table-3.csv')


FC_th<- c(c(0.6,0.8),seq(1,8,0.5)) # log2FC threshold

df <- data.frame(log2FC = numeric(), N_gens = numeric(), N_Intersect_genes = numeric()) # Initilize a data frame store results

for (i in FC_th) {
  
  Atx2_RIP_UP<-DEG_Atx2[which(DEG_Atx2$log2FoldChange > i),'ext_gene'] # RIP IP up regulated genes for each FC_th
  common_genes<-length(intersect(Atx2_RIP_UP,Intersect_gens$ext_gene)) # Intersect genes of up regulated and sleep/wake
  
  # Create a new row as a data frame
  new_row <- data.frame(log2FC = i, N_gens = length(Atx2_RIP_UP), N_Intersect_genes = common_genes) # save all three columns
  
  # Append the new row to the data frame
  df <- rbind(df, new_row)
  
}

# Double axis plot

coeff=0.5 # trick to scale down the yaxis -2 (change the value based on your need)

library(ggplot2)

ggplot(df,aes(x=log2FC))+
  geom_point( aes(y=N_gens),color='red')+
  geom_point(aes(y=N_Intersect_genes/coeff),color='blue')+
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Number of enriched Atx2 RIP IP gense",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(transform =~.*coeff, name="Number of Atx2 RIP genes that intersect with sleep/wake")
  ) +
  theme(
    axis.title.y = element_text(color = 'red', size=10),
    axis.title.y.right = element_text(color = 'blue', size=10)
  )




