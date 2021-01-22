library(goseq)

exonlength_filename = '../data/exon_length.tsv'
edger_filename      = 'output/intersect.tsv'
genes_kept_filename = 'output/keep.intersect.txt'

# load in table listing length of exonic region of each gene

map             = read.table(exonlength_filename, header=TRUE)
rownames(map)   = map$ensembl

# load in edger output

diff = read.table(edger_filename, sep='\t', header=TRUE)
rownames(diff) = diff$ensembl

# load in genes used by edgeR

keep = scan(genes_kept_filename, what="", sep="\n")

# prepare variables

assayed.genes = unique(keep)
de.genes      = rownames(diff)

gene.vector = as.integer(assayed.genes%in%de.genes)
names(gene.vector) = assayed.genes

gene_lengths = map[keep,'exons_length']

# run goseq

pwf=nullp(gene.vector, bias.data = gene_lengths)
GO.wall=goseq(pwf,"hg19","ensGene")

# get FDR for over and under represented terms

GO.wall$over_rep_FDR = p.adjust(GO.wall$over_represented_pvalue, method="BH")
GO.wall$under_rep_FDR = p.adjust(GO.wall$under_represented_pvalue, method="BH")

# keep significantly overrepresented terms

enriched.GO= GO.wall[GO.wall$over_rep_FDR < .05,]

# save table
write.table(enriched.GO, 'output/go_enrichment.tsv', sep='\t', row.names=F)

