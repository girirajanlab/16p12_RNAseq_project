#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=phaser
#SBATCH -o logs/phaser_%a.log
#SBATCH -e logs/phaser_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=4G
#SBATCH --chdir /data5/16p12_RNA/scripts/ase_no_replicates
#SBATCH --array=1-32
#SBATCH --nodelist qingyu

echo `date` started on $HOSTNAME

phaser_dir=/data5/software/phaser
tmp_dir=/data5/16p12_RNA/scripts/ase/tmp


subject=`cat subjects.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`
bam=bams_combined/$subject.sorted.bam
family=`cat subject2vcf.list | grep $subject | cut -f1 | head -n1`
vcf=../ase/vcfs/whatshap.${family}.combined.chrless.vcf.gz

echo `date` $subject
echo `date` $bam
echo `date` $family
echo `date` $vcf

python $phaser_dir/phaser/phaser.py \
	--vcf $vcf \
	--bam $bam \
	--paired_end 1 \
	--mapq 255 \
	--baseq 10 \
	--sample $subject \
	--blacklist hg19_hla.bed \
	--haplo_count_blacklist hg19_haplo_count_blacklist.bed \
	--threads 4 \
	--o output/${subject} \
	--temp_dir $tmp_dir


echo `date` finished


