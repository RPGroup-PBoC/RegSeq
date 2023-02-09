#!/bin/bash

# Name of filtered sequencing file
IN=$1

PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT=${PARENT_PATH##*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}


IN_PATH=$PARENT_PATH"/data/sequencing_data/filtered_sequencing/"$IN"_barcodes.fastq"
mkdir $PARENT_PATH"/data/sequencing_data/bc_by_gc/"

rm -rf $PARENT_PATH"/data/sequencing_data/bc_by_gc/temp/"
mkdir $PARENT_PATH"/data/sequencing_data/bc_by_gc/temp/"


echo "Filtering..."
cat $IN_PATH | awk '/TATTAGGCTTCTCCTCAGCG/' | awk '/TCACTGGCCGTCGTTTTACATGACTGACTGA/' | awk -v o=$PARENT_PATH"/data/sequencing_data/bc_by_gc/temp" 'FNR==1{++f} \
f==1 {a[$2]=$1} \
f==2 {b[$1]=$2} \
f==3 {ind1=substr($0, 0, 4); bc=substr($0,60,20); group=substr($0, 25, 4); if((group in a) && (ind1 in b)){
printf "%s\n", bc >>o"/"b[ind1]"_"a[group]".txt"}}' $PARENT_PATH/data/index_group.txt $PARENT_PATH/data/gc_barcodes.txt -


echo "Counting unique barcodes..."
OUT=$PARENT_PATH"/data/sequencing_data/bc_by_gc"/temp/*.txt
for FILE in $OUT;do
    filename="${FILE##*/}"
    filename="${filename%.*}"
    sort --parallel 20 -T ./ $FILE | uniq -c | sort --parallel 20  -bgr -T ./|  awk -v OFS="\t" '$1=$1' > $PARENT_PATH"/data/sequencing_data/bc_by_gc/"$filename"_expression_counted.csv";
done
