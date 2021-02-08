#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=bam_combine
#SBATCH -o logs/bam_combine_%a.log
#SBATCH -e logs/bam_combine_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/16p12_RNA/scripts/ase_no_replicates
#SBATCH --array=1-32
#SBATCH --nodelist durga

echo `date` started on $HOSTNAME


subject=`cat subjects.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
samples=`cat sample2subject.list | grep $subject | cut -f1`

sample1=`echo $samples | cut -f1 -d' '`
sample2=`echo $samples | cut -f1 -d' '`
sample3=`echo $samples | cut -f1 -d' '`

bam1=/data5/16p12_RNA/scripts/picard/output/${sample1}.Aligned.sortedByCoord.out.md.bam
bam2=/data5/16p12_RNA/scripts/picard/output/${sample2}.Aligned.sortedByCoord.out.md.bam
bam3=/data5/16p12_RNA/scripts/picard/output/${sample3}.Aligned.sortedByCoord.out.md.bam

output=bams_combined/$subject.bam

echo `date` $subject $sample1 $sample2 $sample3
echo `date` $bam1
echo `date` $bam2
echo `date` $bam3
echo `date` $output

# samtools merge bams_combined/$subject.bam $bam1 $bam2 $bam3
samtools sort -m 9G -o bams_combined/$subject.sorted.bam bams_combined/$subject.bam
samtools index bams_combined/$subject.sorted.bam

echo `date` fin