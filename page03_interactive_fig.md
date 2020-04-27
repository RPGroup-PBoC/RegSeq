---
layout: page
title: Interactive Figure
permalink: interactive_a
sidebar: true
interactive: main_interactive_fig.html
---
---

## Figure Description
This website provides interactive figures of information footprints, energy matrices
and regulatory cartoons for all of the genes and growth conditions considered in our experiments.
The gene and growth condition are chosen using the drop down boxes on the left. The
gene regulatory architecture, with identified transcription factors, and their
DNA-protein energy matrix for each transcription factor is displayed below
the information footprint.

The regulatory cartoons are meant to be coarse grained representations of the
regulation. For example, an activator upstream of an RNAP in the regulatory
cartoon indicates an activator upstream of an RNAP site in the wild-type gene,
but the exact distance between the RNAP site and the cartoon is not significant.

The maximum information value for the y-axis displayed in the figure is set by the
maximum value of the 95 percent confidence intervals for the information values
 for any of the displayed base pairs for the gene and growth condition combination.
 If the dataset for the gene and growth condition are unavailable then no
information footprint will be displayed. If a repressor or activator cartoon are
displayed without a protein name, then the identity of the binding site is
unknown.

<!-- The below line includes the interactive figure. Do not change! -->
<center>

{% include_relative interactives/{{page.interactive}} %}

</center>
