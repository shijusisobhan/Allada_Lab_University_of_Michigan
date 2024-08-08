
rm(list=ls())

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/Neuro_transmitters_Receptors')


# Identify neuro transmitters in the BT neurons

# Load Neuro transmitter,receptor and peptide list

NT<-read.csv('NeurotransmitterMarkers.csv')
NR<-read.csv('TMreceptors.csv')
NP<-read.csv('Neuropeptides.csv')

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