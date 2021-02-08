library(peer)

tpm_filename = '16p12.2_rnaseq_analysis/data/16p12_lcl_gene_tpm_by_subject.tsv'
tpm = read.table(tpm_filename, sep='\t', header=T)
rownames(tpm) = tpm$ensembl
tpm$ensembl = NULL

keep_filename = '16p12.2_rnaseq_analysis/outlier_expression_analysis/keep.tsv'
keep = read.table(keep_filename, sep='\t', header=T)

tpm = tpm[rownames(tpm) %in% keep$ensembl,]
tpm = t(tpm)

tpmlog = log2(tpm + 1)

tpmscale = tpmlog

for (i in 1:ncol(tpmscale)) {
	tpmscale[, i] = scale(tpmscale[, i], center=TRUE, scale=TRUE)
}

write.table(tpmscale, 'peer/tpm.scaled.tsv', col.names=F, sep='\t')


model = PEER()
PEER_setPhenoMean(model,as.matrix(tpmscale))
PEER_setNk(model,1)
PEER_update(model)

factors = PEER_getX(model)
resid   = PEER_getResiduals(model)

rownames(resid) = rownames(tpm)
colnames(resid) = colnames(tpm)

rownames(factors) = rownames(tpm)

write.table(factors, 'peer/tpm.factors.1.tsv', col.names=F, sep='\t')
write.table(resid, 'peer/tpm.resid.1.tsv',sep='\t')








