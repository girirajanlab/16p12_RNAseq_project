#/bin/bash


# bed 2 vcf 
/data5/software/plink --recode vcf --snps-only --bfile plink_files/16p12_RNA_genotype --out plink_files/16p12_RNA_vcf

# reformat vcf
jupyter nbconvert --execute  prepare_vcf_by_replicate.ipynb

# vcf 2 bed
/data5/software/plink --vcf 16p12_RNA_by_replicate.vcf --make-bed --out plink_files/16p12_RNA_by_rep

# reformat fam file
jupyter nbconvert --execute prepare_fam_by_replicate.ipynb

# gemma
/data5/software/gemma-0.98.3-linux-static -bfile plink_files/16p12_RNA_by_rep -gk -o relatedness_matrix_by_replicate

