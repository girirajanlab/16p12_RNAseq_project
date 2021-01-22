library(corrplot)

# load in data

df = read.table('../data/all_gene_expression_differences.tsv', header=T)

row.names(df) = df$subject
df$subject = NULL

# get pearson correlations and p-values

corrs = cor(df)
pvals = cor.mtest(df)

# save to figure

pdf('correlation_plot.pdf')

corrplot(corrs, p.mat = pvals$p, 
         insig = "label_sig",  
         pch.col = "white",
        type='upper')

dev.off()
