# Sequence processing scripts

In this folder you can find bash scripts that perform the processing of sequence data and put it in the format needed to perform MCMC step to infer information footprints.
To perform quality filtering on sequence files, we use [`fastp`](https://github.com/OpenGene/fastp). Follow the installation instructions and create a conda environment called `fastp`.

## `1_1_filter_mapping_reads.sh`
Perform quality filtering on the `fastq` sequencing files obtained for identifying promoter variants and associated barcodes. It is assumed that SRA10971993 from the SRA database has been downloaded and split into individual reads, stored in `data/sequencing_data/raw_sequencing/`. Outputs `fastq` files that pass the filters in `data/sequencing_data/filtered_sequencing/`.

## `1_2_extract_promoter_mapping.sh`
Takes filtered sequencing data as input, extracts promoter and barcode pairs from each sequence and counts how often each pair is found in the data. Output is stored by gene group in `data/sequencing_data/mapping_sequences/`.

## `2_1_filter_expression_reads.sh`
Perform quality filtering on the `fastq` sequencing files obtained for counting barcodes in various growth conditions. It is assumed that the respective file has been obtained from the SRA database and stored in `data/sequencing_data/raw_sequencing/`. Outputs `fastq` files that pass the filters in `data/sequencing_data/filtered_sequencing/`.

## `2_2_demultiplex_expression.sh`
Extract barcodes from reads and sort them by index. There are two indices, one gives the growth condition and the other the gene group. Barcodes are sorted into files depending on the indices, sorted and counted. Outputs are stored in `data/sequencing_data/bc_by_gc`.

## `2_3_identify_barcodes.sh`
Takes files containing counted barcodes as input and associates the respective promoter for each found barcode.

