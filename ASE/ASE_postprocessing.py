import pandas as pd
import numpy as np

from statsmodels.stats.multitest import multipletests
from scipy.stats import combine_pvalues
from scipy.stats import binom_test

# load in data from Phaser Gene AE output

ase_file = '../data/ase_gene_ae.tsv'
output_file = 'fisher_combine.tsv'

df = pd.read_csv(ase_file, sep='\t')

# Binomial Test for Allele Specific Expression

a = df.aCount.to_numpy()
b = df.bCount.to_numpy()

p = np.vectorize(binom_test)(a, a+b, p=0.5)

df['p'] = p

# Fisher combine p-values

ref_df = []

df['subject_name'] = df.subject + df.name

for sn in df.subject_name.unique():
    pvalues = df[df.subject_name == sn].p.to_list()
    
    subject = sn[:5]
    gene    = sn[5:]
    
    s, p_combined = combine_pvalues(pvalues)
    
    app = [subject, gene, p_combined]
    ref_df.append(app)
    
ref_df = pd.DataFrame(ref_df, columns = ['subject', 'gene', 'p_combined'])

# for each sample, do multiple testing correction

ref_df['FDR'] = '.'

for sub in ref_df.subject.unique():
    sub_df = ref_df[ref_df.subject == sub].copy()
    fdr = multipletests(sub_df.p_combined)[1]
    ref_df.loc[sub_df.index, 'FDR'] = fdr

# Save file to output

ref_df.to_csv(output_file, sep='\t', index=False)
