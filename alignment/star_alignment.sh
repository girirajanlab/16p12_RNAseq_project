#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=star
#SBATCH -o logs/submit_%a.log
#SBATCH -e errs/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --workdir /data5/16p12_RNA/scripts/star_alignment
#SBATCH --array 1-96%5
#SBATCH --nodelist qingyu

echo `date` starting job on $HOSTNAME

star_index_dir=/data5/16p12_RNA/scripts/star_index/output
fastq_dir=/data5/16p12_RNA/scripts/star_alignment/trim/output2

sample=`cat samples.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`

echo $sample

fastq1="${sample}_R1_trim.fq.gz"
fastq2="${sample}_R2_trim.fq.gz"

python3 run_STAR.py \
        $star_index_dir \
        ${fastq_dir}/$fastq1 \
        ${fastq_dir}/$fastq2 \
        $sample  \
        --threads 4 \
        --output_dir output \
        --hostname $HOSTNAME



echo `date` "done"
