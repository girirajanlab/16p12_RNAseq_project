import pandas as pd
import sys

# arguments
variant_file = sys.argv[1]
network_file = sys.argv[2]
output_file = sys.argv[3]

# python3 rare_var2exp_change.py all_gene_network_simple.tsv /data5/16p12_RNA/scripts/16p12.2_rnaseq_analysis/perm_network/permute_100/random_gene_network_95.tsv output/test.tsv


# load in outlier and variant information
outs = pd.read_csv(variant_file, sep='\t')

# load in network
net = pd.read_csv(network_file, sep='\t')
net.columns = ['startg', 'endg', 'distance',
       'num_connector_genes', 'connector_genes']

# remove genes where start gene is end gene
net = net[net.startg != net.endg]

# get list of all samples
subjects = list(outs['Sample'].unique())

# get list of genes in network
network_genes = list(set(net.startg.to_list()))

# create df to hold output statistics
df = pd.DataFrame(index=subjects)
df['avg_shortest_distance'] = 0.0



# list of gene expression changes
cols_gene_exp_changes = ['Outlier expr.', 'Diff. expression', 'Alt. splicing', 'ASE', 'eQTL allele']
cols_variants         = ['LOF', 'Splicing', 'Missense', 'STR exonic', 'Dup. encapsulated', 'Del. encapsulated']


trios = ['P1C_01','FC_01','P1C_04','P2C_04','SNC_04','P1C_05',
		 'MC_05','P2C_05','P1C_07','SNC_07','M1C_07','P2C_07',
		 'P3C_07','M2C_07','P2C_52','S2NC_52','S3NC_52']

# variants to gene expression changes
for sub in subjects:
	# if subject is not a trio child then skip
	if sub not in trios:
		continue
		
    # print(sub)
    # get all records for sample
	sub_outs = outs[outs['Sample'] == sub]
    
    # get variants
	genes_with_variants = []
	for col in cols_variants:
		append_genes = sub_outs[sub_outs[col] > 0]['Gene name'].to_list()
		genes_with_variants = genes_with_variants + append_genes
		genes_with_variants = list(set(genes_with_variants))
        
    # get gene expression changes
	genes_with_expression_changes = []
	for col in cols_gene_exp_changes:
		append_genes = sub_outs[sub_outs[col] != '.']['Gene name'].to_list()
		genes_with_expression_changes = genes_with_expression_changes + append_genes
		genes_with_expression_changes = list(set(genes_with_expression_changes))
        
    # filter out genes not in network
	genes_with_variants = [s for s in genes_with_variants if s in network_genes]
	genes_with_expression_changes = [s for s in genes_with_expression_changes if s in network_genes]
    
    # get network for variants/expression changes
	subnet = net[net.startg.isin(genes_with_variants)]
	subnet = subnet[subnet.endg.isin(genes_with_expression_changes)]
    
	avg_shortest_distance = subnet.distance.mean()
	df.at[sub, 'avg_shortest_distance'] = avg_shortest_distance


# remove the extra samples from dataframe
df = df[df.avg_shortest_distance != 0]
df.to_csv(output_file, sep='\t', index=True)