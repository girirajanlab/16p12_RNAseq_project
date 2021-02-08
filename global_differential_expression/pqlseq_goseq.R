library(goseq)
suppressPackageStartupMessages(library('GenomicFeatures'))

map    = read.table('../16p12.2_rnaseq_analysis/differential_expression_analysis/exon_length.tsv', header=TRUE)
rownames(map)   = map$ensembl


diff = read.table('pqlseq_by_replicate_FDR.tsv', sep='\t', header=TRUE)
rownames(diff) = diff$ensembl
keep = read.table('pqlseq_by_replicate.tsv', sep='\t', header=TRUE)
keep = rownames(keep)

assayed.genes = unique(keep)
de.genes = rownames(diff)

gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
table(gene.vector)

gene_lengths = map[keep,'exons_length']
head(gene_lengths)

pwf=nullp(gene.vector, bias.data = gene_lengths)

go=goseq(pwf,"hg19","ensGene")

go$over_rep_FDR = p.adjust(go$over_represented_pvalue, method="BH")
go$under_rep_FDR = p.adjust(go$under_represented_pvalue, method="BH")

enriched = go[go$over_rep_FDR < .05,]
dim(enriched)
enriched

filename = 'go_enrichment.tsv'
write.table(enriched, filename, sep='\t', row.names=F)
