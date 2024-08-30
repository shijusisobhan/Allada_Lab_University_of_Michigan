#!/bin/bash


# Define the path to the transcriptome index file
TRANSCRIPTOME_INDEX=/nfs/turbo/umms-rallada/Sequencing_data/transcripts_V_0_50_1.idx

# Define the directory where the fastq files are located
FASTQ_DIR=/nfs/turbo/umms-rallada/Sequencing_data/HY_10969/fastqs_10969-HY

# Define the output directory for Kallisto analysis results
OUTPUT_DIR=/nfs/turbo/umms-rallada/Sequencing_data/HY_10969/fastqs_10969-HY/Kallisto_out

# Loop through each sample and perform Kallisto analysis
samples=(
    "ZT0-A 10969-HY-1"
    "ZT12-A 10969-HY-2"
    "SD-A 10969-HY-3"
    "ZT0-B 10969-HY-4"
    "ZT12-B 10969-HY-5"
    "SD-B 10969-HY-6" 
)


for sample in "${samples[@]}"; do
    sample_name=$(echo $sample | cut -d ' ' -f 1)
    sample_prefix=$(echo $sample | cut -d ' ' -f 2)
    echo "Processing sample: $sample_name"

    # Find fastq files for the current sample
    fastq_files=$(find $FASTQ_DIR -type f -name "${sample_prefix}_*.gz")

    # Run Kallisto analysis
    kallisto quant -i $TRANSCRIPTOME_INDEX -o $OUTPUT_DIR/$sample_name -b 50 $fastq_files
done
