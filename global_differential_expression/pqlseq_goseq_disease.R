library(goseq)
suppressPackageStartupMessages(library('GenomicFeatures'))

map    = read.table('../16p12.2_rnaseq_analysis/differential_expression_analysis/exon_length.tsv', header=TRUE)
rownames(map)   = map$ensembl

diff = read.table('pqlseq_by_replicate_FDR.tsv', sep='\t', header=TRUE)
rownames(diff) = diff$ensembl

anno = read.table('ndd_gene_cateogires.tsv', header=T, sep='\t', stringsAsFactors = F)

split_to_vec = function(s) {
    if (s == ''){
        return(c())
    }
    return(unlist(strsplit(s, ';', fixed=T)))
}
go_map = lapply(anno$database, split_to_vec)
names(go_map) = anno$ensembl


keep = anno[,'ensembl']
diff = diff[diff$ensembl %in% keep, ]

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

# GO.wall=goseq(pwf,"hg19","ensGene")

write.table(go, 'disease_genes_enrichment.tsv', sep='\t', row.names=F)