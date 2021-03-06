# expression quantitative trait loci (eQTL)

### tpm_by_subject.ipynb

Get median TPM for each subject in RNAseq cohort.

### filter_cov.ipynb

Filter genes for at least 50% of samples with > 0.1 TPM.

### format_beds.ipynb

Format bed files for input into QTLtools

### qtltools_pca.sh

Run PCA analysis on genotypes and gene expression of RNAseq cohort.

### combine_covariates.ipynb

Combine explicit (sex, carrier status, family) and principal component covariates.

### qtltools_perm.sh

Run QTLtools on RNAseq cohort.

### qtltools_FDR.sh

Concat output of QTLtools cis output and use QTLtools script runFDR_cis.R to get FDR values and filter for eQTLs with FDR < 0.05.

### annotate_eqtls.ipynb

Annotate eQTLs and filter for eQTLs that are associated with protein coding genes.

### summary.ipynb

Get summary of eQTL minor alleles in RNAseq cohort.

## GTEx Minor eQTLs

Scripts to get number of GTEx minor eQTLs in 16p12 RNA cohort.

### prepare_annovar_input.ipynb

Prepares GTEx eQTL for input into annovar.

### annovar.sh

Runs annovar on GTEx eQTLs to get gnomad frequencies.

### get_gt_info.ipynb

Get GT information for each person in the 16p12 RNA cohort.

### get_maf_eqtl_freq.ipynb

Summarize number of GTEx minor eQTLs per person.

### get_maf_eqtl_freq_disease.ipynb

Summarize number of GTEx minor eQTLs per person by disease gene set category.







