#!/bin/bash

# paths to reads
read1=$1
read2=$2

# Find working directiory
result=${PWD##*/}

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Go back path
parent_path=${parent_path%/*}

# Find data directory
folder=$parent_path'/data/sequencing_data/'

# Output
output=$folder'processed_sequencing/dna_mapping.fastq.gz'

# Process reads
fastp --in1 $read1 --in2 $read2 --merge --merged_out $output --verbose --overlap_len_require 1

# Take out sequences from files
gunzip -c $folder'processed_sequencing/dna_mapping.fastq.gz' | awk 'NR%4==2 {print $1}' > $folder'processed_sequencing/dna_mapping_sequences.txt'