#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=trim
#SBATCH -o logs2/submit_%a.log
#SBATCH -e logs2/submit_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --workdir /data5/16p12_RNA/scripts/star_alignment/trim
#SBATCH --array 1-97

echo `date` starting job on $HOSTNAME

out_dir=/data5/16p12_RNA/scripts/star_alignment/trim/output2
fastq_dir=/data/16p12_LCL_RNASeq/NovaSeq_0519/G459-N059
samples_dir=/data5/16p12_RNA/scripts/star_alignment/sample_fastq

sample=`cat sample.list2 | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`

echo $sample

fastq1="${sample}_R1_001"
fastq2="${sample}_R2_001"

FWD=${fastq_dir}/${fastq1}.fastq.gz
REV=${fastq_dir}/${fastq2}.fastq.gz


java -jar /data/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 $FWD $REV ${out_dir}/${sample}_R1_trim.fq.gz ${out_dir}/${sample}_R1_unpaired.fq.gz ${out_dir}/${sample}_R2_trim.fq.gz ${out_dir}/${sample}_R2_unpaired.fq.gz MINLEN:30




echo `date` "done"
