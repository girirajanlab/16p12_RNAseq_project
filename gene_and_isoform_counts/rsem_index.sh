#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=rsem_index
#SBATCH -o submit.log
#SBATCH -e submit.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --workdir /data5/16p12_RNA/scripts/rsem_index

echo `date` starting job on $HOSTNAME

output=/data5/16p12_RNA/scripts/rsem_index/output
rsem_dir="/data5/anastasia/sw/RSEM-1.2.22"
ref_file=/data/bx_references/GRCh37/human_g1k_v37.fasta
gencode_annotations=/data5/16p12_RNA/GENCODE/gencode.v19.annotation.patched_contigs.gtf


${rsem_dir}/rsem-prepare-reference --num-threads 4 \
	--gtf $gencode_annotations \
	$ref_file \
	output/rsem_reference


echo `date` "done"
