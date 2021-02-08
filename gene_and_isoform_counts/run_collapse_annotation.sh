#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=gene_model
#SBATCH -o submit.log
#SBATCH -e submit.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --workdir /data5/16p12_RNA/scripts/gene_model

echo `date` starting job on $HOSTNAME

ref_file=/data/bx_references/GRCh37/human_g1k_v37.fasta
gencode_annotations=/data5/16p12_RNA/GENCODE/gencode.v19.annotation.patched_contigs.gtf

python3 collapse_annotation.py  $gencode_annotations output/gencode.v19.genes.patched_contigs.gtf --transcript_blacklist gencode19_unannotated_readthrough_blacklist.txt 

# python3 collapse_annotation.py gencode.v26.GRCh38.annotation.gtf gencode.v26.GRCh38.genes.gtf

echo `date` "done"
