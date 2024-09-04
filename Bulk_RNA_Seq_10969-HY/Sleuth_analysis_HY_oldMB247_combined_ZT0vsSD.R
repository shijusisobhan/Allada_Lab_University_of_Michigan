
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

models(so)
