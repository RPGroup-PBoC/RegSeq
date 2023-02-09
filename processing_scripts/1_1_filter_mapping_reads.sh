#!/bin/bash
# Find working directiory
PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}

IN1=$PARENT_PATH"/data/sequencing_data/raw_sequencing/SRR10971993_1.fastq"
IN2=$PARENT_PATH"/data/sequencing_data/raw_sequencing/SRR10971993_2.fastq"
OUT=$PARENT_PATH"/data/sequencing_data/filtered_sequencing/mapping_merged.fastq"

mkdir $PARENT_PATH"/data/sequencing_data/filtered_sequencing"

HTML=$PARENT_PATH"/mapping_fastp_report.html"
JSON=$PARENT_PATH"/mapping_fastp_report.json"

source activate fastp

fastp --in1 $IN1 --in2 $IN2 --merged_out $OUT  --verbose --html $HTML --json $JSON --report_title $html_report --thread '6' --merge --overlap_len_require '3' --n_base_limit '0' 