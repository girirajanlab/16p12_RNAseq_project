# Alignment of RNAseq 

Following the [GTEx Consortium RNAseq pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq).

### Software used

- Trimmomatic 0.36
- STAR 2.4.2a
- Picard 2.9.0
- Python 3.7
- [tin.py 3.0.1](http://rseqc.sourceforge.net/#tin-py)

### trimmomatic.sh

Trims reads less than 30 nucleotide bases long.

### star_index.sh

Generates STAR index of hg19 reference genome.

### star_alingment.sh

Aligns RNAseq to reference genome by running run_STAR.py.

### picard.sh

Marks duplicate reads.

### tin.sh

Calculates transcript integrity scores.
