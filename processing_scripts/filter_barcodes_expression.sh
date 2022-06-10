#!/bin/bash

# index for files
group=$1

# Find working directiory
result=${PWD##*/}

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# Go back path
parent_path=${parent_path%/*}

# Find data directory
folder=$parent_path'/data/processed_sequencing/'

# Output
read=$folder'processed_sequencing/SRR10838518.fastq'
output=$folder'processed_sequencing/SRR10838518_filtered.fastq'

fastp -i $read -o $output --verbose -q 20

cat SRR10838518_filtered.fastq | awk 'NR%4==2 {print $0}' | awk /TGACTGACTG/ > SRR10838518_filtered.txt