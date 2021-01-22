import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# input filenames

edger_filename            = 'output/intersect.tsv'
keep_genes_filename       = 'output/keep.intersect.txt'
singlecell_filename       = '../data/single_cell_expression.tsv'
sc_descriptions_filename  = '../data/single_cell_descriptions.tsv'

# load in data

de       = pd.read_csv(edger_filename, sep='\t')
de_genes = de.ensembl.tolist()

with open(keep_genes_filename, 'r') as f:
    keep = f.readlines()
keep = [s.strip() for s in keep]


desc = pd.read_csv(sc_descriptions_filename, sep='\t')
desc = desc.drop_duplicates('WGCNAcluster')
desc = desc.set_index('WGCNAcluster', drop=False)

te = pd.read_csv(singlecell_filename, sep='\t')
te = te.set_index('ensembl', drop=True)

# keep only genes used by edgeR

te    = te[te.index.to_series().isin(keep)].copy()
te_de = te[te.index.to_series().isin(de_genes)].copy()

total_genes    = te.shape[0]
total_genes_de = te_de.shape[0]

# get counts for number of times DE genes appear in tissues

stats              = pd.DataFrame(index=te.columns)
stats['tissue']    = te.columns
stats['count_all'] = 0
stats['count_de']  = 0

for col in te.columns:
    stats.at[col, 'count_all'] = te[col].sum()
    stats.at[col, 'count_de']  = te_de[col].sum()

# do Fisher's exect test for enrichment

for tissue in stats.index:
    total_genes_not_de = total_genes - total_genes_de
    total_genes_tissue = stats.at[tissue, 'count_all']
    total_genes_not_tissue = total_genes - total_genes_tissue
    
    count_de_tissue = stats.at[tissue, 'count_de']
    count_de_not_tissue = total_genes_de - count_de_tissue
    count_not_de_tissue = total_genes_tissue - count_de_tissue
    count_not_de_not_tissue = total_genes_not_tissue - count_de_not_tissue
    
    p = fisher_exact([[count_de_tissue, count_de_not_tissue], [count_not_de_tissue, count_not_de_not_tissue]])[1]
    stats.at[tissue, 'p'] = p

# multiple testing correction

stats['FDR'] = multipletests(stats['p'], method='fdr_bh')[1]

# add cell type descriptions to table and drop cells with None as clustering info

if 'None' in stats.index:
    stats = stats.drop('None')
stats['description'] = desc.loc[stats.index]['cell type']

# save file

stats.to_csv('output/single_cell_enrichment_of_de_genes.tsv', sep='\t', index=False)