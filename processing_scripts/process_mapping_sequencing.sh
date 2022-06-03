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

# Files
read1=$folder'raw_sequencing/SRR10971993_1.fastq'
read2=$folder'raw_sequencing/SRR10971993_2.fastq'

# Output
output=$folder'processed_sequencing/dna_mapping.fastq.gz'

fastp --in1 $read1 --in2 $read2 --merge --merged_out $output --verbose --overlap_len_require 1

gunzip -c $folder'processed_sequencing/dna_mapping.fastq.gz' | awk 'NR%4==2 {print $1}' > 'processed_sequencing/dna_mapping_sequences.txt'