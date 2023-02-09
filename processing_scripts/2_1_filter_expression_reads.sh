#!/bin/bash
# Find working directiory
PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
SRA=$1
OUT=$2
IN_PATH=$PARENT_PATH"/data/sequencing_data/raw_sequencing/$SRA.fastq"
OUT_PATH=$PARENT_PATH"/data/sequencing_data/filtered_sequencing/"$OUT"_barcodes.fastq"

HTML=$PARENT_PATH"/data/sequencing_data/filtered_sequencing/mapping_fastp_report.html"
JSON=$PARENT_PATH"/data/sequencing_data/filtered_sequencing/mapping_fastp_report.json"

source activate fastp
fastp -i $IN_PATH -o $OUT_PATH  --verbose --html $HTML --json $JSON --report_title $html_report --thread '12'  --average_qual '30'  --n_base_limit '0' --unqualified_percent_limit '10'
