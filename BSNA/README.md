# Brain-Specific Network Analysis

Scripts to get connector genes and shortest distance metrics in a [brain-specific network](https://doi.org/10.1038/nn.4353).

### networkX.ipynb

Get shortest distances from weighted Brain-specific network.

### de_genes2de_genes.py 

Get shortest distance between global DE genes in a network.

### de_perm.sh

Submit de_genes2de_genes.py as a batch job.

### individual_perm.sh

Get shortest distance between an individuals rare variants and gene expression changes in a network.

### rare_var2exp_change.py

Submit rare_var2exp_change.py as a batch job.


# OLD [REMOVE]

### global_de_connectors.ipynb

Get connector genes between all global DE genes.

### global_de_connectors_panther.ipynb

GO term enrichment analysis of global DE connector genes.

### connectors_outlier_and_rare_variant.ipynb

Get connector genes between outliers and rare variants

### connectors_outlier_and_rare_variant_go.ipynb

Run GO term analysis on outlier and rare variant connectors.

### permute_network.ipynb

Create permuted network based on brain-specific network.

### outlier2rare_variant.ipynb

Calculate mean shortest distances between outlier genes and rare variant genes in a brain-specific network.

### outlier2rare_variant_permuted.ipynb

Calculate mean shortest distances between outlier genes and rare variant genes in a permuted network.

### de2rare_variant.ipynb

Calculate mean shortest distances between DE genes and rare variant genes in a brain-specific network.

### de2rare_variant_permuted.ipynb

Calculate mean shortest distances between DE genes and rare variant genes in a permuted network.

### globalde2globalde.ipynb

Calculate mean shortest distances between global DE genes in a brain-specific network.

### globalde2globalde_permuted.ipynb

Calculate mean shortest distances between global DE genes in a permuted network.

### permuted_network_box_plots.ipynb

Create box plots comparing brain-specific network and permuted network (Figures 5E, 5F, S4B)

### conf_intervals_r.R

Get confidence intervals for differences between brain-specific network and permuted network.

