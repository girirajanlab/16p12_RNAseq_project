#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=picard
#SBATCH -o logs/submit_%a.log
#SBATCH -e logs/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --workdir /data5/16p12_RNA/scripts/picard
#SBATCH --array 1-96%5

echo `date` starting job on $HOSTNAME

ref_file=/data/bx_references/GRCh37/human_g1k_v37.fasta
gencode_annotations=/data5/16p12_RNA/GENCODE/gencode.v19.annotation.patched_contigs.gtf

sample=`cat samples.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
star_alignment=/data5/16p12_RNA/scripts/star_alignment/output/${sample}.Aligned.sortedByCoord.out.bam
output_bam=output/${sample}.Aligned.sortedByCoord.out.md.bam
stats_out=${sample}.marked_dup_metrics.txt

java -jar /data/software/picard-tools-2.9.0/picard.jar MarkDuplicates \
	I=$star_alignment \
	O=$output_bam \
	PROGRAM_RECORD_ID=null \
	M=$stats_out \
	ASSUME_SORT_ORDER=coordinate \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100

	

echo `date` "done"
