#!/bin/bash

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
file=$folder'processed_sequencing/dna_mapping.txt'
echo $file
# There are two possible sequence lengths, 299 and 295
# Extract barcode and promoter sequences
# Count unique combinations to reduce file size when identifying promoters
#cat $file | awk '{ if ((length($1) == 295) || (length($1) == 299)) print substr($1, 256, 20)" "substr($1, 21, 160);}' | sort | uniq -c | sort -bgr > $folder'barcode_mapping/dna_mapping_collapsed.csv'
cat $file | awk '{if ((length($0) == 295)) {
    seq=substr($0, 21, length($0)-40);
    print substr(seq, 1, 160)" "substr(seq, length(seq)-19, 20);
        }
}' | sort | uniq -c | sort -bgr > $folder'barcode_mapping/dna_mapping_collapsed.csv'