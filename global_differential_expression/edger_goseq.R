# goseq

library(goseq)
suppressPackageStartupMessages(library('GenomicFeatures'))

map    = read.table('../data/exon_length.tsv', header=TRUE)
rownames(map)   = map$ensembl

diff = read.table('output/no_sex/intersect.tsv', sep='\t', header=TRUE)
keep = scan('output/no_sexcd/keep.intersect.txt', what="", sep="\n")

rownames(diff) = diff$ensembl

assayed.genes = unique(keep)
de.genes = rownames(diff)

gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
table(gene.vector)

gene_lengths = map[keep,'exons_length']
head(gene_lengths)

pwf=nullp(gene.vector, bias.data = gene_lengths)

GO.wall=goseq(pwf,"hg19","ensGene")

print('here')

GO.wall$over_rep_FDR = p.adjust(GO.wall$over_represented_pvalue, method="BH")

enriched.GO= GO.wall[GO.wall$over_rep_FDR < .05,]
dim(enriched.GO)

write.table(GO.wall, 'output/no_sex/intersect.go.tsv', sep='\t', row.names=F)

