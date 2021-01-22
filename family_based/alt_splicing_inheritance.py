#This script identifies inheritance patterns of isoforms with alternative splicing among trios with the 16p12.1 deletion.

#USE: The sample IDs of the prband, carrier parent, and non-carrier parent are input from the command line:
#	python alt_splicing_inheritance.py SGXXX SGYYY SGZZZ
#The script will locate filtered (FDR<0.05) DESeq2 files for each comparison (proband/carrier and proband/noncarrier) in the same directory.

#The output files contain isoforms and logFC values that exhibit alternative splicing (A) de novo, (B) inherited from carrier, and
#(C) inherited from non-carrier.

import sys
import csv

#Get sample IDs from command line
proband=sys.argv[1]
carrier=sys.argv[2]
noncarrier=sys.argv[3]

#Open files for each trio comparison
infile1=open(proband+"_"+carrier+"_DESeq2_filter.csv",'r')
inlines1=infile1.readlines()[1:]
infile2=open(proband+"_"+noncarrier+"_DESeq2_filter.csv",'r')
inlines2=infile2.readlines()[1:]

#Parse gene lists from both comparisons into separate dictionaries
genes_1={} #Changes from carrier parent
genes_2={} #Changes from non-carrier parent

for line1 in inlines1:
	line1=line1.strip().split(",")
	gene1=line1[0].strip('"')
	log2fc1=line1[1]
	genes_1[gene1]=log2fc1

for line2 in inlines2:
	line2=line2.strip().split(",")
	gene2=line2[0].strip('"')
	log2fc2=line2[1]
	genes_2[gene2]=log2fc2

#Identify genes appearing in both lists as "de novo" genes
genes_both={}
for gene1 in genes_1:
	if gene1 in genes_2: #De novo change
		log2fc1=genes_1[gene1]
		log2fc2=genes_2[gene1]
		genes_both[gene1]=[log2fc1,log2fc2]

#Remove de novo genes from initial dictionaries
for denovo_gene in genes_both:
	del genes_1[denovo_gene]
	del genes_2[denovo_gene]

#Output three list files: de novo (appears in both lists), inh. carrier (appears in gene list from non-carrier), 
#and inh. non-carrier (appears in gene list from carrier)

outwriter_denovo = csv.writer(open(proband+"_denovo.txt",'w'),delimiter='\t')
outwriter_denovo.writerow(["ID","log2FC_carrier","log2FC_noncarrier"])
for gene_denovo in genes_both:
	outrow=[gene_denovo,genes_both[gene_denovo][0],genes_both[gene_denovo][1]] #Output both fold-changes
	outwriter_denovo.writerow(outrow)

outwriter_carrier = csv.writer(open(proband+"_inh_carrier.txt",'w'),delimiter='\t')
outwriter_carrier.writerow(["ID","log2FC_carrier","log2FC_noncarrier"])
for gene_carrier in genes_2:
	outrow=[gene_carrier,".",genes_2[gene_carrier]]
	outwriter_carrier.writerow(outrow)

outwriter_noncarrier = csv.writer(open(proband+"_inh_noncarrier.txt",'w'),delimiter='\t')
outwriter_noncarrier.writerow(["ID","log2FC_carrier","log2FC_noncarrier"])
for gene_noncarrier in genes_1:
	outrow=[gene_noncarrier,genes_1[gene_noncarrier],"."]
	outwriter_noncarrier.writerow(outrow)
