# Global differential expression

Code for global differential expression analysis section of paper.


### edgeR.ipynb

This script is used to calculate differential expression between carriers and non-carriers of the 16p12.1 deletion. As input, it takes the raw read counts of the full cohort and the carrier status of all of the samples in the cohort.

This script also runs differential expression between carriers and non-carriers excluding one individual at a time.

### intersection_exclude_one.ipynb

Takes the intersect of gene differentially expressed in our cohort with each sample excluded.

### edger_goseq.R

Runs go enrichment analysis on the differentially expressed genes using goseq.

### edger_goseq_disease.R

Runs go enrichment analysis on the differentially expressed genes in five disease gene categories using goseq.

### edger_goseq_brainspan.R

Runs go enrichment analysis on the differentially expressed genes in brain span preferentially expressed tissues.

### edger_goseq_single_cell.R

Runs go enrichment analysis on the differentially expressed genes in five disease gene categories in genes preferentially expressed in single cell clusters.

### edger_goseq_gtex.R

Runs go enrichment analysis on the differentially expressed genes in GTEx preferentially expressed genes.

### prepare_vcf_by_replicate.ipynb

For GEMMA, prepare vcf file by taking in old vcf file and duplicating sample records x3.

### prepare_fam_by_replicate.ipynb

For GEMMA, prepare fam file by taking in old fam file and renaming samples by replicate names.

### GEMMA.sh

Get relatedness matrix between replicates using WGS data. Runs prepare_vcf_by_replicate.ipynb and prepare_fam_by_replicate.ipynb internally.

### pqlseq.ipynb

Run PQLseq relatedness corrected differential expression between carriers and non-carriers.

### edger_pqlseq_overlap.ipynb

Overlap between edger and PQLseq (Figure S2A).

### pqlseq_goseq.R

Runs go enrichment analysis on the differentially expressed genes found by PQLseq using goseq.

### pqlseq_goseq_disease.R

Runs go enrichment analysis on the differentially expressed genes found by PQLseq in five disease gene categories using goseq.

