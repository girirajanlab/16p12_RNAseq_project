import pandas as pd
import sys

# arguments
de_file = sys.argv[1]
network_file = sys.argv[2]
output_file = sys.argv[3]

# defile=../differential_expression_analysis3/output/no_sex/intersect.tsv
# networkfile=/data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network/permute_100/random_gene_network_${SLURM_ARRAY_TASK_ID}.tsv
# outputfile=output/de_perm_${SLURM_ARRAY_TASK_ID}.tsv
# python3 de_genes2de_genes.py $variantfile $networkfile $outputfile


# load in DE genes
de = pd.read_csv(de_file, sep='\t')
de = de.set_index('ensembl', drop=False)


# load in network
net = pd.read_csv(network_file, sep='\t')
net.columns = ['startg', 'endg', 'distance',
       'num_connector_genes', 'connector_genes']

# remove genes where start gene is end gene
net = net[net.startg != net.endg]

# get list of genes in network
network_genes = list(set(net.startg.to_list()))

# get list of DE genes
de_genes = de.gene.to_list()

# filter out genes not in network
de_genes = [s for s in de_genes if s in network_genes]

# get network for variants/expression changes
subnet = net[net.startg.isin(de_genes)].copy()
subnet = subnet[subnet.endg.isin(de_genes)].copy()

# get average shortest distance and save to file
avg_shortest_distance = subnet.distance.mean()

with open(output_file, 'w') as f:
    f.write(str(avg_shortest_distance))
    f.write('\n')