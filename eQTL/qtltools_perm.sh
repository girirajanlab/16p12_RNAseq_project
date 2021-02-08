#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=perm
#SBATCH -o logs/submit_%a.log
#SBATCH -e logs/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --workdir /data5/16p12_RNA/scripts/eqtl
#SBATCH --array=1-25
#SBATCH --nodelist qingyu


echo `date` started on $HOSTNAME
 
chrom=`cat chrom.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`

echo `date` $chrom
/data5/software/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64 cis --vcf /data3/16p12_WGS/phasing/whatshap/combined.vcf.gz --bed 16p12_lcl_gene_tpm.by_subject.50perc.anno.bed.gz --cov covariates.combined.txt --permute 1000 --out permutations.${chrom}.txt --region $chrom


echo `date` fin


