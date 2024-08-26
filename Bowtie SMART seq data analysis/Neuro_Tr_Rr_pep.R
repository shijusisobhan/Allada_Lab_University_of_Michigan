
rm(list=ls())

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/Neuro_transmitters_Receptors')


# Identify neuro transmitters in the BT neurons

# Load Neuro transmitter,receptor and peptide list

NT<-read.csv('NeurotransmitterMarkers.csv')
NR<-read.csv('TMreceptors.csv')
NP<-read.csv('Neuropeptides.csv')
rm(list=ls())

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/Neuro_transmitters_Receptors')


# Identify neuro transmitters in the BT neurons

# Load Neuro transmitter,receptor and peptide list

NT<-read.csv('NeurotransmitterMarkers.csv')
NR<-read.csv('TMreceptors.csv')
NP<-read.csv('Neuropeptides.csv')

# Load DEG of BT scRNA seq exp

DEG_BT<-read.csv('Normalized_single_cell_count_matrix.csv')
DEG_BT<-DEG_BT[rowSums(DEG_BT[-1] >=1) > 0,]


DEG_MB247<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BTvsMB247_ZT0_ZT12.csv')
MB247_up<-DEG_MB247[which(DEG_MB247$padj<0.05 & DEG_MB247$log2FoldChange>1),'X']

DEG_R5<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_allvsR5_MC.csv')
R5_up<-DEG_R5[which(DEG_R5$padj<0.05 & DEG_R5$log2FoldChange>1),'ext_gene']

DEG_R85<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BTvsR85_ZT0.csv')
R85_up<-DEG_R85[which(DEG_R85$padj<0.05 & DEG_R85$log2FoldChange>1),'X']

DEG_vGAT<-read.csv('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BTvsvGAT_ZT0.csv')
vGAT_up<-DEG_vGAT[which(DEG_vGAT$padj<0.05 & DEG_vGAT$log2FoldChange>1),'X']


DEG_Ma<-read.csv('Ma_RosBash_scRNA_seq_DEG.csv')
DEG_Ma <-DEG_Ma[DEG_Ma$p_val_adj<0.1,]

DEG_Ma_LPN <- subset(DEG_Ma, grepl(c("LPN"), cluster))
DEG_Ma_DN3 <- subset(DEG_Ma, grepl(c("DN3"), cluster))



columns_list<-list(
  
  NT_BT<-DEG_BT[DEG_BT$ext_gene %in% NT$Gene,"ext_gene" ],
  NR_BT<-DEG_BT[DEG_BT$ext_gene %in% NR$Gene,"ext_gene"  ],
  NP_BT<-DEG_BT[DEG_BT$ext_gene %in% NP$Genes, "ext_gene" ],
  
  N.transmitter.LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NT$Gene,"Gene"],
  N.Receptor.LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NR$Gene, "Gene"],
  N.Peptide.LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NP$Gene,"Gene" ],
  
  N.transmitter.DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NT$Gene,"Gene"],
  N.Receptor.DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NR$Gene, "Gene"],
  N.Peptide.DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NP$Gene,"Gene" ],
  
  NT_BTvsMB247<-MB247_up[MB247_up %in% NT$Gene],
  NR_BTvsMB247<-MB247_up[MB247_up %in% NR$Gene],
  NP_BTvsMB247<-MB247_up[MB247_up %in% NP$Gene],
  
  NT_BTvsR5<-R5_up[R5_up %in% NT$Gene],
  NR_BTvsR5<-R5_up[R5_up %in% NR$Gene],
  NP_BTvsR5<-R5_up[R5_up %in% NP$Gene],
  
  NT_BTvsR85<-R85_up[R85_up %in% NT$Gene],
  NR_BTvsR85<-R85_up[R85_up %in% NR$Gene],
  NP_BTvsR85<-R85_up[R85_up %in% NP$Gene],
  
  NT_BTvsvGAT<-vGAT_up[vGAT_up %in% NT$Gene],
  NR_BTvsvGAT<-vGAT_up[vGAT_up %in% NR$Gene],
  NP_BTvsvGAT<-vGAT_up[vGAT_up %in% NP$Gene]
  
  
)


# Find the maximum length of the columns
max_length <- max(sapply(columns_list, length))

# Pad shorter columns with NA to match the maximum length
columns_list <- lapply(columns_list, function(col) {
  c(col, rep(NA, max_length - length(col)))
})

# Combine columns into a data frame
df <- data.frame(columns_list)

colnames(df)<- c('N.transmitter.Bowtie','N.Receptor.Bowtie','N.Peptide.Bowtie',
                 'N.transmitter.DEG.LPN','N.Receptor.DEG.LPN','N.Peptide.DEG.LPN',
                 'N.transmitter.DEG.DN3','N.Receptor.DEG.DN3','N.Peptide.DEG.DN3',
                 'N.transmitter.BTvsMB247.up','N.Receptor.BTvsMB247.up','N.Peptide.BTvsMB247.up',
                 'N.transmitter.BTvsR5.up','N.Receptor.BTvsR5.up','N.Peptide.BTvsR5.up',
                 'N.transmitter.BTvsR85.up','N.Receptor.BTvsR85.up','N.Peptide.BTvsR85.up',
                 'N.transmitter.BTvsvGAT.up','N.Receptor.vGAT.up','N.Peptide.BTvsvGAT.up')

# Write the data frame to a CSV file
write.csv(df, "Neuro_Tr_Rr_Pt_Expressed_in_BT.csv", row.names = FALSE)

# Load DEG of BT scRNA seq exp

DEG_BT<-read.csv('Single_Cell_DEG-BT.csv')
DEG_BT<-DEG_BT[DEG_BT$D_KL>=2,]

#DEG_BT<-DEG_BT[DEG_BT$Pval<0.1,]


NT_BT<-DEG_BT[DEG_BT$Gene %in% NT$Gene, ]
NR_BT<-DEG_BT[DEG_BT$Gene %in% NR$Gene, ]
NP_BT<-DEG_BT[DEG_BT$Gene %in% NP$Gene, ]


# Identify neuro transmitters in the LPN,DN3(Ma & Rosbash)
# Load DEG Rsobash data

DEG_Ma<-read.csv('Ma_RosBash_scRNA_seq_DEG.csv')
DEG_Ma <-DEG_Ma[DEG_Ma$p_val_adj<0.1,]

DEG_Ma_LPN <- subset(DEG_Ma, grepl(c("LPN"), cluster))
DEG_Ma_DN3 <- subset(DEG_Ma, grepl(c("DN3"), cluster))

NT_Ma_LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NT$Gene, ]
NR_Ma_LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NR$Gene, ]
NP_Ma_LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NP$Gene, ]

NT_Ma_DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NT$Gene, ]
NR_Ma_DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NR$Gene, ]
NP_Ma_DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NP$Gene, ]


columns_list<-list(
  
  N.transmitter.Bowtie<-DEG_BT[DEG_BT$Gene %in% NT$Gene, "Gene"],
  N.Receptor.Bowtie<-DEG_BT[DEG_BT$Gene %in% NR$Gene, "Gene"],
  N.Peptide.Bowtie<-DEG_BT[DEG_BT$Gene %in% NP$Gene, "Gene"],
  
  N.transmitter.LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NT$Gene,"Gene"],
  N.Receptor.LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NR$Gene, "Gene"],
  N.Peptide.LPN<-DEG_Ma_LPN[DEG_Ma_LPN$Gene %in% NP$Gene,"Gene" ],
  
  N.transmitter.DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NT$Gene,"Gene"],
  N.Receptor.DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NR$Gene, "Gene"],
  N.Peptide.DN3<-DEG_Ma_DN3[DEG_Ma_DN3$Gene %in% NP$Gene,"Gene" ]
  
)

# Find the maximum length of the columns
max_length <- max(sapply(columns_list, length))

# Pad shorter columns with NA to match the maximum length
columns_list <- lapply(columns_list, function(col) {
  c(col, rep(NA, max_length - length(col)))
})

# Combine columns into a data frame
df <- data.frame(columns_list)

colnames(df)<- c('N.transmitter.Bowtie','N.Receptor.Bowtie',
                 'N.Peptide.Bowtie', 'N.transmitter.LPN',
                 'N.Receptor.LPN', 'N.Peptide.LPN',
                 'N.transmitter.DN3', 'N.Receptor.DN3', 'N.Peptide.DN3')

# Write the data frame to a CSV file
write.csv(df, "Neuro_Transmitter_Receptor_peptide.csv", row.names = FALSE)