#!/bin/bash
gc=$1
group=$2


PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}

# File containing barcode-promoter key
IN_MAP=$PARENT_PATH"/data/sequencing_data/mapping_sequences/"${group}"_mapping_identified.csv"
# File containing barcode counts
IN_DNA_BC="/home/tom/git/1000_genes_ecoli/data/bc_by_gc/"${gc}"_DNA_"${group}"_expression_counted.csv"
IN_RNA_BC="/home/tom/git/1000_genes_ecoli/data/bc_by_gc/"${gc}"_RNA_"${group}"_expression_counted.csv"
# location to store identified barcodes
OUT_DNA="/home/tom/git/1000_genes_ecoli/data/bc_by_gc/"${gc}"_DNA_"${group}"_identified.txt"
OUT_RNA="/home/tom/git/1000_genes_ecoli/data/bc_by_gc/"${gc}"_RNA_"${group}"_identified.txt"

awk '((NR==FNR) && (NR>1)){a[$1]=$0; next} ($2 in a){b=$2;$2="";print $0" "a[b]}' $IN_MAP $IN_DNA_BC  > $OUT_DNA
awk '((NR==FNR) && (NR>1)){a[$1]=$0; next} ($2 in a){b=$2;$2="";print $0" "a[b]}' $IN_MAP $IN_RNA_BC  > $OUT_RNA


#join -j 2 -a 1 -a 2 -e 0 -o 1.1,2.1,1.2,1.3,1.4,1.5,1.6 temp_RNA.txt temp_DNA.txt > temp.txt


#echo -e "ct_1 ct_0 barcode promoter mapping_count name nmut "| cat - "temp.txt"  > $OUT_BC

#rm "temp.txt"
#rm "temp_DNA.txt"
#rm "temp_RNA.txt"