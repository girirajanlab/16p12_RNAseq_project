#!/bin/bash

# run annovar on GTEx eQTLs


perl /data4/software/annovar/table_annovar.pl gtex_eqtls_annovar_input.tsv  /data4/software/annovar/humandb/ \
 -buildver hg19 \
 -out gtex_eqtls \
 -remove \
 -protocol gnomad_genome \
 -operation f \
 -nastring . \
 -thread 6 \
 -xref /data4/software/annovar/humandb/gene_annotations.txt 
 
 
 