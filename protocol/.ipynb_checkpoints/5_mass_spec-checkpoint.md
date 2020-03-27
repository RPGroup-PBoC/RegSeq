# RegSeq Mass Spectrometry

## DNA Affinity Chromatography


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

Uncertainty due to number of independent sequences

1400 promoter variants were ordered from TWIST Bioscience for each promoter studied. Due
to errors in synthesis, additional mutations are introduced into the ordered oligos. As a result,
the final number of variants received was an average of 2200 per promoter. To test whether
the number of promoter variants is a significant source of uncertainty in the experiment we
computationally reduced the number of promoter variants used in the analysis of the zapAB -10
RNAP region. Each sub-sampling was performed 3 times. The results, as displayed in Figure S3,
show that there is only a small effect on the resulting sequence logo until the library has been
reduced to approximately 500 promoter variants.

In some cases, we used an alternative approach to mass spectrometry to discover the TF identity
regulating a given promoter based on sequence analysis using a motif comparison tool. TOMTOM
[14] is a tool that uses a statistical method to infer if a putative motif resembles any previously
discovered motif in a database. Of interest, it accounts for all possible offsets between the motifs.
Moreover, it uses a suite of metrics to compare between motifs such as Kullback-Leibler divergence, Pearson correlation, euclidean distance, among others

We performed comparisons of the motifs generated from our energy matrices to those generated from all known transcription factor binding sites in RegulonDB. Figure S4 shows a result of
TOMTOM, where we compared the motif derived from the -35 region of the ybjX promoter and
found a good match with the motif of PhoP from RegulonDB.
The information derived from this approach was then used to guide some of the TF knockout
experiments, in order to validate its interaction with a target promoter characterized by the loss of
the information footprint. Furthermore, we also used TOMTOM to search for similarities between
our own database of motifs, in order to generate regulatory hypotheses in tandem. This was
particularly useful when looking at the group of GlpR binding sites found in this experiment.