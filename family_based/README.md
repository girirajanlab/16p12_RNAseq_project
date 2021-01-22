# Family Based Analysis

Scripts used for the family based analysis of the paper.

### Libraries used

R 3.6.3:

* [edgeR 3.28.1](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* tximport
* DESeq2

Python 3.7.3:

* sys
* csv


### edger.R

This script is used to calculate differential expression between parent-child pairs. As input, it takes the raw read counts of the full cohort and a file listing the family relationship between samples.

### DESeq2_tximport_isoform.r

This script calculates differential expression between parent-child pairs using DESeq2 on an isoform level.

### DESeq2_tximpot_gene.r

This script calculates differential expression between parent-child pairs using DESeq2 on a gene level.

### isoform_to_gene.py

Converts ensembl isoform names to gene names.

### alt_splicing_inheritance.py

Script that identifies inheritance patterns for genes with alternative splicing.
