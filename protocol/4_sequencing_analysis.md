# RegSeq Sequence Design

## Initial Tidying of Sequence Files and Counting of Barcode Frequencies

## Information Footprints
We use information footprints as a tool for hypothesis generation to identify regions which may
contain transcription factor binding sites. In general, a mutation within a transcription factor site
is likely to severely weaken that site. We look for groups of positions where mutation away from
wild type has a large effect on gene expression. Our data sets consist of nucleotide sequences, the
number of times we sequenced the construct in the plasmid library, and the number of times we
sequenced its corresponding mRNA. A simplified data set on a 4 nucleotide sequence then might
look like

| Sequence    | DNA Counts  | mRNA Counts |
| ----------- | ----------- | ----------- |
| ACTA        | 5           | 23          |
| ATTA        | 5           | 3           |
| CCTG        | 11          | 11          |
| TAGA        | 12          | 3           |
| GTGC        | 2           | 0           |
| CACA        | 8           | 7           |
| AGGC        | 7           | 3           |

One possible calculation to measure the impact of a given mutation on expression is to take
all sequences which have base b at position i and determine the number of mRNAs produced per
read in the sequencing library. By comparing the values for different bases we could determine
how large of an effect mutation has on gene expression. However, in this paper we will use mutual
information to quantify the effect of mutation, as [4] demonstrated could be done successfully.
In Table 1 the frequency of the different nucleotides in the library at position 2 is 40% A, 32% C, 14% G and 14% T. Cytosine is enriched in the mRNA transcripts over the original library, as it
now composes 68% of all mRNA sequencing reads while A, G, and T only compose only 20%,
6%, and 6% respectively. Large enrichment of some bases over others occurs when base identity
is important for gene expression. We can quantify how important using the mutual information
between base identity and gene expression level. Mutual information is given at position i by



## Insertion of Barcoded Library into Plasmid Backbone

## Cloning or Integration of Barcoded Library

## Pooling of _E. coli_ and DNA / RNA Isolation


To determine putative transcription factor binding sites, we first compute the effect of mutations
on gene expression at a base pair-by-base pair level using information footprints. The information
footprints are a hypothesis generating tool and we choose which regions to further investigate
using techniques such as mass spectrometry by visually inspecting the data for regions of 10
to 20 base pairs that have high information content compared to background. Our technique
currently relies on using human intuition to determine binding sites, but to validate these choices
and to capture all regions important for gene expression we computationally identify regions
where gene expression is changed significantly up or down by mutation (p < 0.01), and discard
any potential sites which do not fit this criteria. We infer the effect of mutation using Markov
Chain Monte Carlo, and we use the distribution of parameters from the inference to form a 99 %
confidence interval for the average effect of mutation across a 15 base pair region. We include
binding sites that are statistically significant at the 0.01 level in any of the tested growth conditions.
Many false positives will be secondary RNAP sites and we remove from consideration any
sites that resemble RNAP sites. We fit energy matrices to each of the possible binding sites and
use the preferred DNA sequence for binding to identify the RNAP sites. We use both visual
inspection to compare the preferred sequence to known consensus sequences for each of the E.
coli sigma factor binding sites (for example, do the preferred bases in the energy matrix have
few mismatches to the TGNTATAAT extended minus 10 for σ
70 sites), and the TOMTOM tool
[47] to computationally compare the potential site to examples of σ
70
, σ
38, and σ
54 sites that we
determined in this experiment. For further details see Supplementary Figure S4. We discard any
sites that have a p-value of similarity with an RNAP site of less than 5x10−3
in the TOMTOM
analysis or are deemed to be too visually similar to RNAP sites. If a single site contains an RNAP
site along with a transcription factor site we remove only those bases containing the probable
24
RNAP site. This results in 95 identified transcription factor binding regions.
For primary RNAP sites, we include a list of probable sigma factor identities as Supplementary
Table 2. Sites are judged by visual similarity to consensus binding sites. Those sites where the true
sigma factor is unclear due to overlapping binding sites are omitted. Overlapping binding sites
(from multiple TFs or RNAP sites) in general can pose issues for this method. In many cases, looking at growth conditions where only one of the relevant transcription factors is present or active is
an effective way to establish site boundaries and infer correct energy matrices. For sites where no
adequate growth condition can be found, or when a TF overlaps with an RNAP site, the energy
matrix will not be reflective of the true DNA-protein interaction energies. If the TFs in overlapping
sites are composed of one activator and one repressor, then we use the point at which the effect
of mutation shifts from activator-like to repressor-like as a demarcation point between binding
sites. We see a case of a potentially overlooked repressor due to overlapping sites in Figure 4(B),
as there are several repressor like bases overlapping the RNAP -10 site and the effect weakens in
low oxygen growth. However, due to the effect of the RNAP site, when averaged over a potential 15 base pair region, the repressor-like bases do not have a significant effect on gene expression