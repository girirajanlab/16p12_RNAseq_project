#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=index_star
#SBATCH -o submit.log
#SBATCH -e submit.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --workdir /data5/16p12_RNA/scripts/star_index


echo `date` starting job on $HOSTNAME

output_dir=/data5/16p12_RNA/scripts/star_index/output
star_dir="/data5/anastasia/sw/STAR-STAR_2.4.2a/qingyu"
# ref_file=/data/bx_references/hg19/ucsc.hg19.fasta
ref_file=/data/bx_references/GRCh37/human_g1k_v37.fasta
gencode_annotations=/data5/16p12_RNA/GENCODE/gencode.v19.annotation.patched_contigs.gtf


${star_dir}/STAR \
        --runMode genomeGenerate \
        --genomeDir $output_dir \
        --genomeFastaFiles $ref_file \
        --sjdbGTFfile $gencode_annotations \
        --sjdbOverhang 52 \
        --runThreadN 50


echo `date` "done"
