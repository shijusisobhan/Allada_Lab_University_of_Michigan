
# sluth analysis Raw data


rm(list=ls())

library(sleuth)
library("gridExtra")
library("cowplot")
library(biomaRt)
library(devtools)


# change the directory

setwd('X:/Sequencing_data/HY_10969/Data_analysis_shiju/HY_oldMB247_combined_ZT0vsSD')

base_dir <- 'X:/Sequencing_data/HY_10969/Data_analysis_shiju/HY_oldMB247_combined_ZT0vsSD/Kallisto_out'


# Sample name
sample_id <- dir(file.path(base_dir))

# Location of alignment output for each sample
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))

# Read the sample-condtion information table (We have to made it based on experiment)
s2c <- read.table("Sample_condition.txt", header = TRUE, stringsAsFactors=FALSE)

# Add a column that locate the directoies of each sample and condtion. This column must be labeled path
# The user should check whether or not the order is correct. 

s2c <- dplyr::mutate(s2c, path = kal_dirs)
# collect gene names  from Ensembl

mart <- useEnsembl(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl")

# Add genename into the sleuth table
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)



#Quality check


Quality_filter <- function(row, min_reads = 5, min_prop = 0.1)
{
  mean(row >= min_reads) >= min_prop
}


so<- sleuth_prep(s2c,target_mapping = t2g, filter_fun=Quality_filter, aggregation_column = 'ext_gene', gene_mode = TRUE, extra_bootstrap_summary=T,
                 read_bootstrap_tpm = TRUE)


# The analysis
# we will fit two models.
# 1. reduced model -includes the batch
# 2. full model - Include 'Batch' and 'condition'
# compare the full model to the reduced model with the likelihood ratio test.
# Identifies genes whose abundances are significantly better explained when 'condition' is taken into account

so <- sleuth_fit(so, ~Batch, 'reduced')
so <- sleuth_fit(so, ~Batch+condition, 'full')
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

plot_pca(so, color_by = 'Batch', text_labels = F, units = "tpm")

plot_pca(so, color_by = 'condition', text_labels = F, units = "tpm")


# 1.How do we know what models are available for testing?
models(so)

beta='conditionSD'

so <- sleuth_wt(so, which_beta = beta)

sig_level=0.05
sig_th =0.6
test_table <- sleuth_results(so, beta)
test_table$raw_b <- exp(test_table$b)
test_table$log2_b <- log2(test_table$raw_b)
test_table$neg_log10_qval<- -log10(test_table$qval)


test_table$diffexpressed <- "Not sig"
# if log2Foldchange > 0 and qvalue < 0.1, set as "UP" 
test_table$diffexpressed[test_table$log2_b > sig_th  & test_table$qval < sig_level] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
test_table$diffexpressed[test_table$log2_b < -sig_th  & test_table$qval < sig_level] <- "DOWN"


N_significant<-length(test_table$diffexpressed[test_table$diffexpressed !="Not sig"])
N_UP<-length(test_table$diffexpressed[test_table$diffexpressed =="UP"])
N_DOWN<-length(test_table$diffexpressed[test_table$diffexpressed =="DOWN"])

N_significant
N_UP
N_DOWN




top_genes<- test_table[1:5, "target_id"]

library('ggplot2')
library(ggrepel)

ggplot(test_table) + geom_point(aes(x = log2_b, y = neg_log10_qval, col=diffexpressed))+
  geom_vline(xintercept=0, col="black",  linetype="dashed")+
  scale_color_manual(values=c("green","black", "red"))+
  labs(x = "log2(FC)", y="-log10(q)", colour="DEG") +
  geom_text_repel(data = subset(test_table, target_id == top_genes), 
                  aes(x = log2_b, y = neg_log10_qval, label = target_id), color = "black")




# get the gene expression table

sleuth_gene_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')

DGE_table<-cbind(row.names(sleuth_gene_matrix), sleuth_gene_matrix)# Add row names as one colum (gene names

p_val_sort<-(test_table[order(test_table$target_id),]) # sort according to gene name


p_val_reduce<-p_val_sort[c("target_id","pval","qval", "log2_b")] # extract only genename,p and q and log2()


DGE_all<-merge(DGE_table, p_val_reduce, by.x=1 , by.y="target_id") # megrge with gene expression and p, q values

DGE_all<-na.omit(DGE_all)


write.csv(DGE_all, 'DEG_Combined_10969-HY_oldMB247_ZT0vsSD.csv')

