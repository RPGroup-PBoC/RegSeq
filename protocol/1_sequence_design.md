# RegSeq Sequence Design

## Finding Regulatory Sites in the _E. coli_ Genome
Genes in this study were chosen to cover several different categories. 29 genes had some information on their regulation already known to validate our method under a number of conditions. 37 were chosen because the work of [1] demonstrated that gene expression changed significantly under different growth conditions. A handful of genes such as minC, maoP, or fdhE were chosen because we found either their physiological significance interesting, as in the case of the cell division gene minC or that we found the gene regulatory question interesting, such for the intraoperon regulation demonstrated by fdhE. The remainder of the genes were chosen because they had no regulatory information, often had minimal information about the function of the gene, and had an annotated transcription start site (TSS) in RegulonDB.

A known limitation of the experiment is that the mutational window is limited to 160 bp. As such, it is important to correctly target the mutation window to the location around the most active TSS. To do this we first prioritized those TSS which have been extensively experimentally validated and catalogued in RegulonDB. Secondly we selected those sites which had evidence of active transcription from RACE experiments [2] and were listed in RegulonDB. If the intergenic region was small enough, we covered the entire region with our mutation window. If none of these options were available, we used computationally predicted start sites.

## Computational Mutations of Regulatory Binding Sites

## Validation of Mutated Regulatory Binding Sites

## Ordering Mutated Oligo Library from TWIST or IDT
Promoter variants were synthesized on a microarray (TWIST Bioscience, San Francisco, CA). The sequences were designed computationally such that each base in the 160 bp promoter region has a 10% probability of being mutated. For each given promoter’s library, we ensured that the mutation rate as averaged across all sequences was kept between 9.5% and 10.5%, otherwise the library was regenerated. There are an average of 2200 unique promoter sequences per gene (for an analysis of how our results depend upon number of unique promoter sequences see Supplementary Figure S3 in the Reg-Seq paper from 2020). An average of 5 unique 20 base pair barcodes per variant promoter was used for the purpose of counting transcripts. 

## Preparation of Oligo Library
Once the oligo library has arrived from Twist (typically the tube arrives with a yellow screw-top lid), use the [Twist DNA Resuspension](https://www.twistbioscience.com/resources/twist-dna-resuspension-guidelines) guidelines to prepare the oligo library. We typically store oligo libraries at a 10 ng / microliter concentration, in 5 microliter aliquots. 




