#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=de_perm
#SBATCH -o logs/de_perm_%a.log
#SBATCH -e logs/de_perm_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network2
#SBATCH --array 1-100
#SBATCH --nodelist ramona

echo `date` starting job on $HOSTNAME

defile=../differential_expression_analysis3/output/no_sex/intersect.tsv
networkfile=/data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network/permute_100/random_gene_network_${SLURM_ARRAY_TASK_ID}.tsv
# networkfile=/data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network/all_gene_network.tsv
outputfile=output/de_perm_${SLURM_ARRAY_TASK_ID}.tsv
# outputfile=output/de_perm_BSNA.tsv

echo $defile
echo $networkfile
echo $outputfile

python3 de_genes2de_genes.py $defile $networkfile $outputfile


echo `date` "done"
