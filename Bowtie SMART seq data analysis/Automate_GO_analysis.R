
library(clusterProfiler)
library(org.Dm.eg.db)  # Use appropriate annotation package for your species 9it is for drosophila melanogaster)
library(pathview)


setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis')

DEG<-read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_allvsWB_ZT0.csv")

gene_list<-DEG[which(DEG$log2FoldChange>0.6 & DEG$padj<0.05), 'ext_gene']
#gene_list_kegg<-paste0('Dmel_',gene_list)

 # df$Name <- paste("Dr.", df$Name)
# Example gene list
#gene_list <- c("Gene1", "Gene2", "Gene3")  # Replace with your gene list

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)

# Perform GO enrichment analysis
ego_BP <- enrichGO(gene         = gene_entrez$ENTREZID,
                OrgDb        = org.Dm.eg.db,
                keyType      = 'ENTREZID',
                ont          = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# Save the results

GO_result_BP<-as.data.frame(ego_BP)
GO_result_BP<-GO_result_BP[,c('Description', 'qvalue', 'geneID', 'Count')]

colnames(GO_result_BP)[colnames(GO_result_BP) == "Description"] <- "Biological_process"




# Perform GO enrichment analysis
ego_CC <- enrichGO(gene         = gene_entrez$ENTREZID,
                   OrgDb        = org.Dm.eg.db,
                   keyType      = 'ENTREZID',
                   ont          = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Save the results

GO_result_CC<-as.data.frame(ego_CC)
GO_result_CC<-GO_result_CC[,c('Description', 'qvalue', 'geneID', 'Count')]

colnames(GO_result_CC)[colnames(GO_result_CC) == "Description"] <- "Cellular_component"




# Perform GO enrichment analysis
ego_MF <- enrichGO(gene         = gene_entrez$ENTREZID,
                   OrgDb        = org.Dm.eg.db,
                   keyType      = 'ENTREZID',
                   ont          = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Save the results

GO_result_MF<-as.data.frame(ego_MF)
GO_result_MF<-GO_result_MF[,c('Description', 'qvalue', 'geneID', 'Count')]

colnames(GO_result_MF)[colnames(GO_result_MF) == "Description"] <- "Molecular_Function"



#*******************************************************************************************
library(gprofiler2)
# Perform enrichment analysis
gost_results <- gost(query = gene_list, 
                     organism = "dmelanogaster",  # organism code for Drosophila melanogaster
                     sources = c("KEGG"))

Kegg_pathway<-data.frame(gost_results$result)

Kegg_pathway<-Kegg_pathway[,]


#*****************************************************************************************

library(ReactomePA)

# Perform Reactome pathway analysis
reactome_results <- enrichPathway(gene = gene_entrez$ENTREZID, 
                                  organism = "fly",  # organism code for Drosophila
                                  pvalueCutoff = 0.05)


Kegg_pathway_results<-data.frame(reactome_results)


library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

addWorksheet(wb, "BP")
writeData(wb, "BP", GO_result_BP)

addWorksheet(wb, "CC")
writeData(wb, "CC", GO_result_CC)

addWorksheet(wb, "MF")
writeData(wb, "MF", GO_result_MF)


# Save the workbook
saveWorkbook(wb, "GO_BT_allvsWB_Up.xlsx", overwrite = TRUE)

