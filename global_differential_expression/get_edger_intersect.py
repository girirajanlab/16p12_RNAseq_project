import pandas as pd

# load in sample list

pheno = pd.read_csv('../data/phenotypes.tsv', sep='\t')
pheno = pheno.drop_duplicates('subject')
pheno = pheno.set_index('subject', drop=False)

# load in edger output for full cohort

orig = pd.read_csv('output/edger_full_cohort.tsv', sep='\t')

# get intersection of edger output with each sample excluded taking note of logFC direction

sub = 'SG001'
df = pd.read_csv('output/edger_{}_excluded.tsv'.format(sub), sep='\t')
df['dir'] = df.logFC.apply(lambda x: '+'if x > 0 else '-')
df['gene_dir'] = df.ensembl + df.dir
intersect = df.gene_dir.to_list()


for sub in pheno.subject:
    df = pd.read_csv('output/edger_{}_excluded.tsv'.format(sub), sep='\t')
    df['dir'] = df.logFC.apply(lambda x: '+'if x > 0 else '-')
    df['gene_dir'] = df.ensembl + df.dir
    excl_de = df.gene_dir.to_list()
    intersect = list(set(intersect) & set(excl_de))
    
intersect = [s[:-1] for s in intersect]

# for intersected genes, use FDR and logFC from edgeR run on full cohort

df_intersect = orig[orig.ensembl.isin(intersect)]

# save intersect table

df_intersect.to_csv('output/intersect.tsv', sep='\t', index=False)

# get intersection of genes used by edger

keep = pd.read_csv('output/keep.SG001.txt', sep='\t', header=None)
keep = keep[0].to_list()
for sub in pheno.subject:
    df = pd.read_csv('output/keep.{}.txt'.format(sub), sep='\t', header=None)
    excl_keep = df[0].to_list()
    keep = list(set(keep) & set(excl_keep))

# save genes used by edgeR in all comparisons
    
with open('output/keep.intersect.txt', 'w') as f:
    for gene in keep:
        f.write(gene)
        f.write('\n')