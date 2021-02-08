#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=pca
#SBATCH -o logs/pca.log
#SBATCH -e logs/pca.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --workdir /data5/16p12_RNA/scripts/eqtl
#SBATCH --nodelist qingyu


echo `date` started on $HOSTNAME
 

/data5/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64 pca --bed 16p12_lcl_gene_tpm.by_subject.50perc.anno.bed.gz --scale --center --out 16p12_lcl_gene_tpm.by_subject.50perc.anno

/data5/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64 pca --vcf /data3/16p12_WGS/phasing/whatshap/combined.vcf.gz --scale --center --maf 0.05 --distance 50000 --out genotypes 



echo `date` fin


