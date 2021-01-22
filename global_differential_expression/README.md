# Global differential expression

Code for global differential expression analysis section of paper.

### Libraries used

R 3.6.3:

* [edgeR 3.28.1](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* [goseq 1.38.0](https://bioconductor.org/packages/release/bioc/html/goseq.html)
* glue 1.4.1

Python 3.7.3:

* pandas 1.0.0
* scipy 1.1.0
* statsmodels 0.11.0

### edger.R

This script is used to calculate differential expression between carriers and non-carriers of the 16p12.1 deletion. As input, it takes the raw read counts of the full cohort and the carrier status of all of the samples in the cohort.

### get_intersect.py

Takes the intersect of gene differentially expressed in our cohort with each sample excluded.

### goseq.R

Runs go enrichment analysis on the differentially expressed genes using goseq.

### GTEx_tissue_enrichment.py

Runs enrichment analysis for genes preferrentially expressed in [GTEx](https://www.gtexportal.org/home/index.html) tissues using a Fisher's exact test.


### BrainSpan_tissue_enrichment.py

Runs enrichment analysis for genes preferrentially expressed in [BrainSpan](https://www.brainspan.org/) tissues at different developmental timepoints using a Fisher's exact test.

### single_cell_tissue_enrichment.py

Runs enrichment analysis for genes preferrentially expressed in [single-cell RNAseq clusters](https://doi.org/10.1126/science.aap8809) using a Fisher's exact test.
