library(goseq)
suppressPackageStartupMessages(library('GenomicFeatures'))

map    = read.table('../data/exon_length.tsv', header=TRUE)
rownames(map)   = map$ensembl

diff = read.table('output/no_sex/intersect.tsv', sep='\t', header=TRUE)
keep = scan('output/no_sex/keep.intersect.txt', what="", sep="\n")

anno = read.table('../tissue_expression2/brainspan_preferential_tissue_enrichment.tsv', header=T, sep='\t', stringsAsFactors = F)

split_to_vec = function(s) {
    if (s == ''){
        return(c())
    }
    return(unlist(strsplit(s, ';', fixed=T)))
}
go_map = lapply(anno$enriched_tissues, split_to_vec)
names(go_map) = anno$ensembl

rownames(diff) = diff$ensembl

# diff = diff[diff$ensembl %in% keep, ]

assayed.genes = unique(keep)
de.genes = rownames(diff)

gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
table(gene.vector)

gene_lengths = map[keep,'exons_length']
head(gene_lengths)

head(gene.vector)



pwf=nullp(gene.vector, bias.data = gene_lengths)

go=goseq(pwf,"hg19","ensGene",gene2cat=go_map, use_genes_without_cat=TRUE)
go$over_rep_FDR = p.adjust(go$over_represented_pvalue, method="BH")


write.table(go, 'output/no_sex/intersect.go.brainspan.tsv', sep='\t', row.names=F)