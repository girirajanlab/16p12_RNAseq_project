#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=genecount
#SBATCH -o logs2/submit_%a.log
#SBATCH -e logs2/submit_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --workdir /data5/16p12_RNA/scripts/gene_counts
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user awt5304@psu.edu
#SBATCH --array=11-96%4
#96


echo `date` starting job on $HOSTNAME

sample=`cat samples.list | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id'`

echo $sample

ref_file=/data/bx_references/GRCh37/human_g1k_v37.fasta
gencode_annotations=../gene_model/output/gencode.v19.genes.patched_contigs.gtf
star_alignment=/data5/16p12_RNA/scripts/picard/output/${sample}.Aligned.sortedByCoord.out.md.bam

mkdir output2/${sample}
output_dir=output2/${sample}/out

/data5/software/jre1.7.0_80/bin/java -jar /data5/software/RNA-SeQC_v1.1.8.jar -strictMode \
	-t $gencode_annotations \
	-r $ref_file \
	-o $output_dir \
	-s $sample,$star_alignment,$sample

	

echo `date` "done"
