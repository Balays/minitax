#!/bin/bash

# Check if the correct number of arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_tsv_file> <destination_directory>"
    exit 1
fi

# Assign arguments to variables
TSV_FILE=$1
DEST_DIR=$2

# Ensure the destination directory exists
mkdir -p "$DEST_DIR"

# Set the number of accessions per batch
BATCH_SIZE=50  # You can adjust this to optimize based on your needs

# Read the TSV file, extract the first column, and process in batches
cut -f1 "$TSV_FILE" | split -l $BATCH_SIZE - "${DEST_DIR}/batch_"

# Loop through each batch and download sequences
for batch_file in "${DEST_DIR}"/batch_*; do
    # Create a comma-separated list of accession numbers for each batch
    accessions=$(tr '\n' ',' < "$batch_file" | sed 's/,$//')

    # Download the batch of sequences in FASTA format
    efetch -db nucleotide -id "$accessions" -format fasta >> "${DEST_DIR}/all_sequences.fasta"

    # Optional: Sleep to avoid rate limiting
    sleep 1
done

# Clean up the temporary batch files
rm "${DEST_DIR}"/batch_*

echo "Download completed."
