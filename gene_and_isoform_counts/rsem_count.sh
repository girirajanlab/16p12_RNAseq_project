#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=rsem
#SBATCH -o logs/submit_%a.log
#SBATCH -e logs/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --workdir /data5/16p12_RNA/scripts/rsem_count
#SBATCH --array 1-96%5

echo `date` starting job on $HOSTNAME

output=/data5/16p12_RNA/scripts/rsem_index/output
rsem_dir="/data5/anastasia/sw/RSEM-1.2.22"
rsem_reference=/data5/16p12_RNA/scripts/rsem_index/output


sample=`cat samples.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
star_alignment=/data5/16p12_RNA/scripts/star_alignment/output/${sample}.Aligned.toTranscriptome.out.bam

python3 run_RSEM.py $rsem_reference \
       $star_alignment \
       output/$sample \
       --threads 4 \
       --hostname $HOSTNAME



echo `date` "done"
