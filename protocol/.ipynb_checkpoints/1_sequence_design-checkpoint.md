# RegSeq Sequence Design

## Finding Regulatory Sites in the _E. coli_ Genome using EcoCyc
The identification of regulatory binding sites -- the first step in the RegSeq protocol -- begins on [EcoCyc](https://ecocyc.org/). Only a limited number of searches can be performed without logging in, so we recommend that you access the [Caltech Library](https://www.library.caltech.edu/) website, select 'Databases' on the right-hand side of the screen, and [search for BioCyc](https://libguides.caltech.edu/az.php?q=biocyc). On the next webpage, log in with your Caltech username and password. The BioCyc website will load. On the home page, select _E. coli_, which will take you to the EcoCyc website.

From the EcoCyc website, you can search for any gene that might be interesting to study. I will outline this process using a sample gene: _marA_. Simply type the gene name, such as marA, in the search box and hit 'Enter'. On the page that loads, there will be several results: Pathways, Genes, Proteins, and so forth. Click on 'marA' under the 'Genes' heading.

This loads a webpage that is devoted to the gene of interest (marA). It provides a lot of useful information, including the gene length, its position on the genome, evidence for the gene's existence (experiment is best), and so forth. The Summary window, in the pane below, shows even more useful information. In particular, it shows that the marA gene is part of an operon (with marR and marB), is transcribed by RNAP σ70, and is regulated by numerous activators (in green boxes) and numerous repressors (in red boxes). Any of the boxes in the figure can be selected, and it will bring you to a page for its gene.

The other useful page to explore, on a gene's page, is the "Operons" tab. This loads a new image which shows the orientations of genes, as well as the milieu of activators and repressors which regulate its expression. In this image, pay attention to a couple of things:

- Any transcription factor (TF) outlined in a SOLID box has better _evidence_ than a TF outlined in a STRIPED/DOTTED box. The latter TFs are often inferred. Hover over any TF to see what experimental evidence, if any, exists for that TF regulating the gene of interest.

- Transcription start sites are given by the bent arrows. Hover over an arrow to see what evidence exists for its position -- a solid arrow, again, has more evidence than a dotted arrow. Hovering over an arrow also provides crucial information on the Transcription Start Site (TSS). If multiple TSS' are shown for a given gene or operon, select the TSS for which the most experimental evidence exists. If these putative TSS' are spread out by ~200bp or more, greater independent judgment must be used to determine which seems reasonable. The best evidence for a TSS position is given by RACE experiments, full datasets for which can be found on [RegulonDB](http://regulondb.ccg.unam.mx/menu/download/datasets/index.jsp).

Carry out your EcoCyc search for all genes of interest. For each gene, record what TFs (if any are known) regulate its expression. Also record the TSS for each gene's operon. These numbers will be important in the next step: generating the mutated libraries.

To get a better understanding for why some genes are chosen in lieu of others, note that, in the original RegSeq paper (2020, Ireland _et al_.), genes were chosen to cover several different categories. 29 genes had some information on their regulation already known to validate the RegSeq method under a number of conditions. 37 genes were chosen because of [previous work](https://www.nature.com/articles/nbt.3418), which demonstrated that gene expression changed significantly under different growth conditions. The remainder of the genes were chosen because they had no regulatory information, often had minimal information about the function of the gene, and had an annotated transcription start site (TSS) in RegulonDB.

A known limitation of the RegSeq experiment is that the mutational window is limited to 160 bp. As such, it is important to correctly target the mutation window to the location around the most active TSS, hence why it is important to choose those TSS that have experimental evidence. Again, active transcription can be determined from [RACE experiments](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007526). If the intergenic region is small enough, the 160bp mutational window may cover the entire area of interest. If none of these options were available (in other words, no experimental evidence exists at all), use computationally predicted start sites.

## Computational Mutations of Regulatory Binding Sites

## Validation of Mutated Regulatory Binding Sites

## Ordering Mutated Oligo Library from TWIST or IDT
Promoter variants were synthesized on a microarray (TWIST Bioscience, San Francisco, CA). The sequences were designed computationally such that each base in the 160 bp promoter region has a 10% probability of being mutated. For each given promoter’s library, we ensured that the mutation rate as averaged across all sequences was kept between 9.5% and 10.5%, otherwise the library was regenerated. There are an average of 2200 unique promoter sequences per gene (for an analysis of how our results depend upon number of unique promoter sequences see Supplementary Figure S3 in the Reg-Seq paper from 2020). An average of 5 unique 20 base pair barcodes per variant promoter was used for the purpose of counting transcripts. 

## Preparation of Oligo Library
Once the oligo library has arrived from Twist (typically the tube arrives with a yellow screw-top lid), use the [Twist DNA Resuspension](https://www.twistbioscience.com/resources/twist-dna-resuspension-guidelines) guidelines to prepare the oligo library. We typically store oligo libraries at a 10 ng / microliter concentration, in 5 microliter aliquots. 




