rm(list=ls())

DEG_Atx2<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Atx2/Writeup/Table-2.csv')

# Wake-sleep genes
# Intersect_gens<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Atx2/Writeup/Atx2_Wake_sleep/Table-3.csv')
Wake_genes<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Atx2/Writeup/Allada_or_Swinderen_Lab_wake_genes_log2b_0.5.csv')

sleep_genes<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Atx2/Writeup/Allada_or_Swinderen_Lab_sleep_genes_log2b_0.5.csv')


# selected genes

Wake_genes<-Wake_genes[which(Wake_genes$frequency>1),'ext_gene']

sleep_genes<-sleep_genes[which(sleep_genes$frequency>1),'ext_gene']

Intersect_gens<-c(Wake_genes,sleep_genes)



FC_th<- c(c(0.6,0.8),seq(1,8,0.5)) # log2FC threshold

df <- data.frame(log2FC = numeric(), N_gens = numeric(), N_Intersect_genes = numeric()) # Initilize a data frame store results

for (i in FC_th) {
  
  Atx2_RIP_UP<-DEG_Atx2[which(DEG_Atx2$log2FoldChange > i),'ext_gene'] # RIP IP up regulated genes for each FC_th
  common_genes<-length(intersect(Atx2_RIP_UP,Intersect_gens)) # Intersect genes of up regulated and sleep/wake
  
  # Create a new row as a data frame
  new_row <- data.frame(log2FC = i, N_gens = length(Atx2_RIP_UP), N_Intersect_genes = common_genes) # save all three columns
  
  # Append the new row to the data frame
  df <- rbind(df, new_row)
  
}

# Double axis plot

coeff=1/6 # trick to scale down the yaxis -2 (change the value based on your need)

library(ggplot2)

ggplot(df,aes(x=log2FC))+
  geom_point( aes(y=N_gens),color='red')+
  geom_point(aes(y=N_Intersect_genes/coeff),color='blue')+
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "# enriched Atx2 RIP IP genes",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(transform =~.*coeff, name="# Atx2 RIP IP enriched genes Intersect with sleep/wake genes")
  ) +
  theme(
    axis.title.y = element_text(color = 'red', size=10),
    axis.title.y.right = element_text(color = 'blue', size=10)
  )




