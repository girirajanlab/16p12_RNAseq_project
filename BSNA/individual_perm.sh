#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=indv_perm
#SBATCH -o logs/indv_perm_%a.log
#SBATCH -e logs/indv_perm_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network2
#SBATCH --array 70-88
#SBATCH --nodelist ramona,durga

echo `date` starting job on $HOSTNAME

variantfile=expression_rare_variant_list.txt
networkfile=/data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network/permute_100/random_gene_network_${SLURM_ARRAY_TASK_ID}.tsv
# networkfile=/data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network/all_gene_network.tsv
outputfile=output/indv_perm_${SLURM_ARRAY_TASK_ID}.tsv
# outputfile=output/indv_perm_BSNA.tsv

echo $variantfile
echo $networkfile
echo $outputfile

python3 rare_var2exp_change.py $variantfile $networkfile $outputfile
# python3 rare_var2exp_change.py expression_rare_variant_list.txt /data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network/permute_100/random_gene_network_${SLURM_ARRAY_TASK_ID}.tsv output/indv_perm_${SLURM_ARRAY_TASK_ID}.tsv


echo `date` "done"
