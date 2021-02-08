#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=tin
#SBATCH -o logs/submit_%a.log
#SBATCH -e logs/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --workdir /data5/16p12_RNA/scripts/rseqc/transcript_integrity
#SBATCH --array 1-96%5


echo `date` starting job on $HOSTNAME

sample=$SLURM_ARRAY_TASK_ID

sample=`head -n $sample ../samples.list | tail -n1`
echo `date` $sample

coord=toTranscriptome
coord=sortedByCoord

bamfile=/data5/16p12_RNA/scripts/star_alignment/output/${sample}.Aligned.${coord}.out.bam

/usr/local/bin/tin.py --input=$bamfile --refgene=hg19_GencodeCompV19.bed > output/${sample}.${coord}.txt
echo `date` finished
