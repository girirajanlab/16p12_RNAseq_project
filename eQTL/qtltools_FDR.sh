#!/bin/bash

##
## Runs QTLtools script runFDR_cis.R 
## https://qtltools.github.io/qtltools/pages/mode_cis_permutation.html
##

# Concat QTLtools cis output
cat ../eqtl/permutations.chr*.txt > permutations.all.txt

# Run runFDR_cis.R
Rscript /data5/software/QTLtools_1.2_Ubuntu16.04_x86_64/script/runFDR_cis.R permutations.all.txt 0.05 permutations.FDR.txt
