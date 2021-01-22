library(edgeR)

# input filenames

relations_filename  = '../data/family_summaries.tsv'
pheno_filename      = '../data/phenotypes.tsv'
gene_reads_filename = '../data/gene_reads.tsv'


pheno           = read.table(pheno_filename, sep='\t', header=TRUE, stringsAsFactors = F)
rownames(pheno) = pheno$sample

relations           = read.table(relations_filename, sep='\t', header=TRUE, stringsAsFactors = F)
rownames(relations) = relations$subject

rawdf               = read.table(gene_reads_filename, sep='\t', header=TRUE)
rownames(rawdf)     = rawdf$Name
rawdf$Name          = NULL

## format data

gencode2ensembl = function(s) {
    return(strsplit(s, '.', fixed=T)[[1]][1])
}

rownames(rawdf) = unlist(lapply(rownames(rawdf), gencode2ensembl))
rawmat = as.matrix(rawdf)

group         = pheno$subject
y             = DGEList(counts=rawmat, group=group)
y             = calcNormFactors(y)
design        = model.matrix(~0+group)

y             = estimateDisp(y, design, robust=TRUE)
fit           = glmQLFit(y,design, robust=TRUE)

subject_map = function(design, subject) {
    subject = paste0('group', subject)
    pos = (1:length(colnames(design)))[colnames(design) == subject]
    return(pos)
}

delta = function(d) {
    if (d < 0) {
        return('-')
    }
    if (d > 0){
        return('+')
    }
}

for (subject1 in rownames(relations)) {

    if (subject1 == 'SG011') next
    for (parent in c('mother', 'father')) {
        subject2 = as.character(relations[subject1,parent])
        print(paste0(subject1, ' ', subject2))
        
        contrast = numeric(length(colnames(design)))
        contrast[subject_map(design, subject1)] = 1
        contrast[subject_map(design, subject2)] = -1
        
        qlf = glmQLFTest(fit,contrast=contrast)
        cdf = topTags(qlf, n=56202, p.value=0.05)$table
        print(dim(cdf))
        
        cdf$direction = as.character(lapply(cdf$logFC, delta))
        
        diff_save = cbind(rownames(cdf), cdf)
        colnames(diff_save)[1] = "ensembl"
        outfile = paste0('output/', subject1, '_', subject2, '.tsv')
        write.table(diff_save, outfile, sep='\t', row.names=F, col.names=T)
        }
    }