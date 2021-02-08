# Gene and Isoform Counts

Following the [GTEx Consortium RNAseq pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq)

### Software used

- RSEM 1.2.22
- RNASeQC 1.1.8
- Python 3.7


### rsem_index.sh

Create index for RSEM based on GENCODE 19 genes.

### rsem_count.sh

Create isoform-level counts.

### run_collapse_annotation.sh

Collapse GENCODE 19 annotation based on genes.

### rnaseqc_count.sh

Create gene-level counts.
