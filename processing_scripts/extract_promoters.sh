#!/bin/bash

# index for files
group=$1

# Find working directiory
result=${PWD##*/}

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Go back path
parent_path=${parent_path%/*}

# Find data directory
folder=$parent_path'/data/sequencing_data/'

# Go to data
cd $folder

# Sequence file names
file='processed_sequencing/dna_mapping_sequences.txt'

# There are two possible sequence lengths, 299 and 295
# Extract barcode and promoter sequences
# Count unique combinations to reduce file size when identifying promoters
cat $file | awk '{ if ((length($1) == 295) || (length($1) == 299)) print substr($1, 256, 20)"\t"substr($1, 21, 160);}' | sort | uniq -c | sort -bgr > 'barcode_mapping/dna_mapping_collape.csv'