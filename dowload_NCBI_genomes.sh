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

# Function to download sequences
download_sequence() {
    accession=$1
    dest_dir=$2
    # Use efetch to download the sequence in FASTA format and save it to the specified directory
    efetch -db nucleotide -id "$accession" -format fasta > "${dest_dir}/${accession}.fasta"
}

export -f download_sequence

# Read the first column from the TSV file, pass each accession to download_sequence function using parallel
cat "$TSV_FILE" | cut -f1 | parallel download_sequence {} "$DEST_DIR"

echo "Download completed."
