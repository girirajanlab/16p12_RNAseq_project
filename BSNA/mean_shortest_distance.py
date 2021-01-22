import pandas as pd


network_filename = '../data/brain_specific_network.tsv'
global_de_filename = '../data/global_de_genes.tsv'

# load input files

net = pd.read_csv(network_filename, sep='\t')
net = net.set_index(['Start gene', 'End gene'], drop=False)

de  = pd.read_csv(global_de_filename, sep='\t')


# function gets mean shortest distance between input gene and a list of genes

def get_mean_shortest_distance(gene, gene_list):
    if gene not in net['Start gene']:
        return 'gene not in network'

    subnet = net.loc[gene].copy()
    subnet = subnet[subnet['End gene'] != gene].copy()
    subnet = subnet[subnet['End gene'].isin(gene_list)].copy()

    mean_distance = subnet['Number of connector genes'].mean()
    return mean_distance

# Find the mean shortest distances between globally differentially expressed genes and all genes in the network

de_genes = de.gene.to_list()
gene_list = list(net['Start gene'].unique()) + list(net['End gene'].unique())
gene_list = list(set(gene_list))

mean_distances = []

for gene in de_genes:
    mean_distance = get_mean_shortest_distance(gene, gene_list)
    mean_distances.append(mean_distance)

# write to file

with open('mean_shortest_distances.txt', 'w') as f:
    for i in range(len(de_genes)):
        gene = de_genes[i]
        mean_distance = mean_distances[i]
        f.write('{}\t{}\n'.format(gene, mean_distance))

