import pandas as pd
import numpy as np

global_de_filename = '../data/global_de_genes.tsv'
network_filename   = '../data/brain_specific_network.tsv'

# load in input files

de  = pd.read_csv(global_de_filename, sep='\t')
net = pd.read_csv(network_filename, sep='\t')

# function that returns connector genes for input gene list

def get_connector_genes(genes):
    
    all_connectors = []
    
    for gene1 in genes:
        subnet = net[net['Start gene'] == gene1]
        subnet = subnet.set_index('End gene', drop=False)
        end_genes = subnet['End gene'].to_list()

        for gene2 in genes:
            if gene2 not in end_genes:
                continue
            connectors = subnet.loc[gene2, 'Connector genes']

            if (type(connectors) == float) or (type(connectors) == np.float64):
                continue

            all_connectors = all_connectors + connectors.split(',')

        all_connectors = list(set(all_connectors))
        
    return all_connectors

# get connector genes for global DE genes

genes = de.gene.to_list()
connectors = get_connector_genes(genes)

# save to file

with open('global_de_connector_genes.txt', 'w') as f:
    for gene in connectors:
        f.write(gene)
        f.write('\n')

