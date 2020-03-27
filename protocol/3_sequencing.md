# RegSeq Sequence Design

## Initial Amplification of Oligo Library
https://www.twistbioscience.com/op_protocol_ampguidelines

## Second Amplification and Barcoding of Oligo Library

## Insertion of Barcoded Library into Plasmid Backbone

## Cloning or Integration of Barcoded Library

## Pooling of _E. coli_ and DNA / RNA Isolation

qPCR was preformed to check the level of DNA contamination and the mRNA tags
were PCR amplified and Illumina sequenced. Within a single growth condition, all promoter
variants for all regulatory regions were tested in a single multiplexed RNA-Seq experiment. All
sequencing was carried out by either the Millard and Muriel Jacobs Genetics and Genomics
Laboratory at Caltech (HiSeq 2500) on a 100 bp single read flow cell or using the sequencing
services from NGX Bio on a 250 bp or 150 base paired end flow cell.


All sequencing was carried out by either the Millard and Muriel Jacobs Genetics and Genomics
Laboratory at Caltech (HiSeq 2500) on a 100 bp single read flow cell or using the sequencing
services from NGX Bio on a 250 bp or 150 base paired end flow cell. The total library was
first sequenced by PCR amplifying the region containing the variant promoters as well as the
corresponding barcodes. This allowed us to uniquely associate each random 20 bp barcode with
a promoter variant. Any barcode which was associated with a promoter variant with insertions
or deletions was removed from further analysis. Similarly, any barcode that was associated with
multiple promoter variants was also removed from the analysis. The paired end reads from this
sequencing step were then assembled using the FLASH tool [3]. Any sequence with PHRED
score less than 20 was removed using the FastX toolkit. Additionally, when sequencing the initial
library, sequences which only appear in the dataset once were not included in further analysis in
order to remove possible sequencing errors.
For all the MPRA experiments, only the region containing the random 20 bp barcode was
sequenced, since the barcode can be matched to a specific promoter variant using the initial library
sequencing run described above. For a given growth condition, each promoter yielded 50,000
to 500,000 usable sequencing reads. Under some growth conditions, genes were not analyzed
further if they did not have at least 50,000 reads.
To determine which base pair regions were statistically significant a 99% confidence interval
was constructed using the MCMC inference to determine the uncertainty.