# Effect size

Scripts perform log odds ratio tests for the association between local rare variants and gene expression changes of three types: outlier genes expression, differential expression between parent and child, and alternative isoforms present in the child but not the parent.

### filter_protein_coding_lcl_expressed.ipynb

Filters genes for protein coding and for at least one subject that has all three replicates with > 0.2 TPM.

### tpm_by_subject.ipynb

Get TPM by subject by taking median TPM of replicates.

### tpm_zscores.ipynb

Get zscores of TPM by subject of protein coding and LCL expressed genes.

### peer.R 

Ruun PEER correction on TPM z-scores. These are the final values we use in the paper.

### outlier_log_odds.py

Calculates log odds ratios for the presence of rare variants near genes with outlier expression.

### diff_exp_log_odds.py

Calculates log odds ratios for the presence of rare variants near genes with differential expression between parent-child pairs.

### alt_splicing_log_odds.py

Calculates log odds ratios for the presence of rare variants near genes with unique isoforms between parent-child pairs.
