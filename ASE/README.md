# Allele Specific Expression

Postproceesing of output from [Phaser Gene AE](https://github.com/secastel/phaser/tree/master/phaser_gene_ae). The script takes gene allelic counts as input and outputs fisher combined and multiple hypothesis corrected p-values for allele specific expression (ASE).

### combine_bams.sh

Combine RNAseq BAM files

### phaser.sh

Run phaser on RNAseq cohort.

### phaser_gene_ae.sh

Run phaser Gene AE on RNAseq cohort.

### binomial_test_and_filter.ipynb

Do binomial test, FDR correction, and filter for protein coding genes.

### cadd_annotation.ipynb

Annotate ASE for deletrious variants.

### summary.ipynb

Summarize ASE by individual

### summary_cadd.ipynb

Summarize ASE with deleterious variant on overexpressed haplotype by individual.
