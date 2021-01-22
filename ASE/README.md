# Allele Specific Expression

Postproceesing of output from [Phaser Gene AE](https://github.com/secastel/phaser/tree/master/phaser_gene_ae). The script takes gene allelic counts as input and outputs fisher combined and multiple hypothesis corrected p-values for allele specific expression (ASE).

### Libraries used

Python 3.7.3:

* pandas 1.0.0
* numpy 1.16.2
* statsmodels 0.11.0
* scipy 1.1.0

### ASE_postprocessing.py

This script takes output from Phaser Gene AE and calculates signficance for ASE for each genes using a binomial test. The script also combines p-values for multiple replicates and does multiple hypothesis correction.
