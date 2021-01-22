#This Python script converts a list of GENCODE hg19 isoforms from DESeq2 analysis to their corresponding genes.
#The input is a text file of isoforms and log-fold changes in each comparison, and a gene:isoform file from GENCODE. 
#The output is the same input text file with gene names instead of isoforms.

import sys
import csv

#Get file name from command line
infile_name=sys.argv[1]

#Open files
infile=open("../"+infile_name+".txt",'r')
inlines=infile.readlines()[1:]
gencode_file=open("gencode_v19_ensembl_NBCI.txt",'r')
gencode_lines=gencode_file.readlines()

#Create dictionary for gene conversion
gene_conversion={}
for gencode_line in gencode_lines:
	gencode_line=gencode_line.strip().split("\t")
	transcript=gencode_line[1]
	gene=gencode_line[0]
	gene_conversion[transcript]=gene

#Convert transcript in input file, and output line to new file
outwriter = csv.writer(open(infile_name+"_genes.txt",'w'),delimiter='\t')
header=["Gene","log2FC"]
outwriter.writerow(header)
for line in inlines:
	line=line.strip().split("\t")
	isoform=line[0]
	log2fc=line[1:]

	new_gene=gene_conversion[isoform]
	outrow=[new_gene]+log2fc
	outwriter.writerow(outrow)