# RegSeq Sequence Design

## Initial Amplification of Oligo Library
https://www.twistbioscience.com/op_protocol_ampguidelines

## Second Amplification and Barcoding of Oligo Library

## Insertion of Barcoded Library into Plasmid Backbone

## Cloning or Integration of Barcoded Library

## Pooling of _E. coli_ and DNA / RNA Isolation

Upon identifying a putative transcription factor binding site, we used DNA affinity chromatography, as done in [20] to isolate and enrich for the transcription factor of interest. In brief, we order
biotinylated oligos of our binding site of interest (Integrated DNA Technologies, Coralville, IA)
along with a control, ”scrambled” sequence, that we expect to have no specificity for the given
transcription factor. We tether these oligos to magnetic streptavidin beads (Dynabeads MyOne T1;
ThermoFisher, Waltham, MA), and incubate them overnight with whole cell lysate grown in the
presences of either heavy (with 15N) or light (with 14N) lysine for the experimental and control
sequences, respectively. The next day, proteins are recovered by digesting the DNA with the PtsI
restriction enzyme (New England Biolabs, Ipswich, MA), whose cut site was incorporated into all
designed oligos.
Protein samples were then prepared for mass spectrometry by either in-gel or in-solution
digestion using the Lys-C protease (Wako Chemicals, Osaka, Japan). Liquid chromatography
coupled mass spectrometry (LC-MS) was performed as previously described by [20], and is
further discussed in the SI. SILAC labeling was performed by growing cells (∆ LysA) in either
heavy isotope form of lysine or its natural form.
It is also important to note that while we relied on the SILAC method to identify the TF identity
for each promoter, our approach doesnt require this specific technique. Specifically, our method
only requires a way to contrast between the copy number of proteins bound to a target promoter
in relation to a scrambled version of the promoter. In principle, one could use multiplexed
proteomics based on isobaric mass tags [55] to characterize up to 10 promoters in parallel. Isobaric
tags are reagents used to covalently modify peptides by using the heavy-isotope distribution in
the tag to encode different conditions. The most widely adopted methods for isobaric tagging are
the isobaric tag for relative and absolute quantitation (iTRAQ) and the tandem mass tag (TMT).
This multiplexed approach involves the fragmentation of peptide ions by colliding with an inert
gas. The resulting ions are resolved in a second MS-MS scan (MS2).
25
Only a subset (13) of all transcription factor targets were identified by mass spectrometry
due to limitations in scaling the technique to large numbers of targets. The transcription factors
identified by this method are enriched more than any other DNA binding protein, with p <
0.01 using the outlier detection method as outlined by Cox and Mann [56], with corrections for
multiple hypothesis testing using the method proposed by Benjamini and Hochberg [57].
