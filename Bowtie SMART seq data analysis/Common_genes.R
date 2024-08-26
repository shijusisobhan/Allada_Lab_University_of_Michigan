

rm(list=ls())
setwd("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis")

DEG_MB247<-read.csv('DEG_BTvsMB247_ZT0_ZT12.csv')
MB247_up<-DEG_MB247[which(DEG_MB247$padj<0.05 & DEG_MB247$log2FoldChange>1),'X']

DEG_R5<-read.csv('DEG_BT_allvsR5_MC.csv')
R5_up<-DEG_R5[which(DEG_R5$padj<0.05 & DEG_R5$log2FoldChange>1),'ext_gene']

DEG_R85<-read.csv('DEG_BTvsR85_ZT0.csv')
R85_up<-DEG_R85[which(DEG_R85$padj<0.05 & DEG_R85$log2FoldChange>1),'X']

DEG_vGAT<-read.csv('DEG_BTvsvGAT_ZT0.csv')
vGAT_up<-DEG_vGAT[which(DEG_vGAT$padj<0.05 & DEG_vGAT$log2FoldChange>1),'X']


library(VennDiagram)
library(RColorBrewer)

# Chart
venn.diagram(
  x = list(MB247_up,R5_up,R85_up,vGAT_up),
  category.names = c('MB247','R5','R85','vGAT'),
  filename = 'venn_UP_247_R5_vgat_R85.png',
  output=TRUE,
  fill = c("green", "red","blue","black"),
)


MB247_R5_R85_VGAT <- Reduce(intersect, list(MB247_up,R5_up,R85_up,vGAT_up))
