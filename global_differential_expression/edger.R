library(edgeR)
library(glue)

## load in data

gene_reads_filename     = '../data/gene_reads.tsv'
pheno_filename          = '../data/phenotypes.tsv'

pheno                   = read.table(pheno_filename, sep='\t', header=TRUE)
rownames(pheno)         = pheno$sample

rawdf               = read.table(gene_reads_filename, sep='\t', header=TRUE)
rownames(rawdf)     = rawdf$Name
rawdf$Name          = NULL

## format data

gencode2ensembl = function(s) {
    return(strsplit(s, '.', fixed=T)[[1]][1])
}

rownames(rawdf) = unlist(lapply(rownames(rawdf), gencode2ensembl))
rawmat = as.matrix(rawdf)
pheno  = pheno[colnames(rawmat),]

## function contrast map creates a vector with group of each sample (for edgeR)

contrast_map = function(design, group) {
    group = paste0('group', group)
    pos = (1:length(colnames(design)))[colnames(design) == group]
    return(pos)
}

## run edger on full cohort

samples = pheno$sample
group   = pheno$status
family  = pheno$family

y = DGEList(counts=rawmat, group=group)
keep = filterByExpr(y)
write.table(names(keep[keep]), glue('output/keep.txt'), 
        row.names = F, col.names = F, quote = F)
y = y[keep,,keep.lib.sizes=FALSE]
y = calcNormFactors(y)

design        = model.matrix(~0+group+family)
rownames(design) = colnames(y)
y = estimateDisp(y, design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)

contrast = numeric(length(colnames(design)))
contrast[contrast_map(design, 'non_carrier')] = -1
contrast[contrast_map(design, 'carrier')] = 1

qlf = glmQLFTest(fit,contrast=contrast)
topdf = topTags(qlf, n=56202, p.value=0.05)$table

## save edger output

outfile = 'output/edger_full_cohort.tsv'
if (is.null(topdf)) {
    write.table(topdf, outfile, sep='\t', row.names=F, col.names=T)

} else{
    diff_save = cbind(rownames(topdf), topdf)
    colnames(diff_save)[1] = "ensembl"
    write.table(diff_save, outfile, sep='\t', row.names=F, col.names=T)
}


## for each sample, exclude that sample and run edgeR on remaining cohort

for (excl_sub in unique(pheno$subject)) {
    
	# remove sample from pheno and gene reads df
    excl_pheno = pheno[pheno$subject != excl_sub,]
    samples = excl_pheno$sample
    excl_rawmat = rawmat[,samples]
    
	# get sample attributes
    group         = excl_pheno$status
    family        = excl_pheno$family
	
	# run edger
    y = DGEList(counts=excl_rawmat, group=group)
    keep = filterByExpr(y)
    write.table(names(keep[keep]), glue('output/keep.{excl_sub}.txt'), 
            row.names = F, col.names = F, quote = F)
    y = y[keep,,keep.lib.sizes=FALSE]
    y = calcNormFactors(y)
    	
    design        = model.matrix(~0+group+family)
    rownames(design) = colnames(y)
    y = estimateDisp(y, design, robust=TRUE)
    fit = glmQLFit(y, design, robust=TRUE)
    
	contrast = numeric(length(colnames(design)))
    contrast[contrast_map(design, 'non_carrier')] = -1
    contrast[contrast_map(design, 'carrier')] = 1
	
    qlf = glmQLFTest(fit,contrast=contrast)
    topdf = topTags(qlf, n=56202, p.value=0.05)$table
    
    # save edgeR output
	outfile = glue('output/edger_{excl_sub}_excluded.tsv')
	if (is.null(topdf)) {
		write.table(topdf, outfile, sep='\t', row.names=F, col.names=T)
	} else{
		diff_save = cbind(rownames(topdf), topdf)
		colnames(diff_save)[1] = "ensembl"
		write.table(diff_save, outfile, sep='\t', row.names=F, col.names=T)
	}

#     diff_save = cbind(rownames(topdf), topdf)
#     colnames(diff_save)[1] = "ensembl"
#     write.table(diff_save, outfile, sep='\t', row.names=F, col.names=T)
	
}