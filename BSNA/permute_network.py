import pandas as pd
from random import choice 

network_filename = '../data/brain_specific_network.tsv'

# load input files

net = pd.read_csv(network_filename, sep='\t')
net = net.set_index(['Start gene', 'End gene'], drop=False)



# make a dict of new gene name mappings

all_genes = list(net['Start gene'].unique())
all_genes = all_genes + list(net['End gene'].unique())
all_genes = list(set(all_genes))

gene2random = {}

genes = all_genes
for gene in genes:
    random = choice(genes)
    gene2random[gene] = random
    
    genes = [s for s in genes if random != s]


# create permuted network

def get_random_gene(gene):
    random = gene2random[gene]
    return random

random_net = net.copy()


random_net['Start gene'] = random_net['Start gene'].apply(get_random_gene)
random_net['End gene'] = random_net['End gene'].apply(get_random_gene)

# write to file

random_net.to_csv('permuted_net.tsv', sep='\t', index=False)

