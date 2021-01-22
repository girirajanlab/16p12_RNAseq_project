# Effect size

Scripts perform log odds ratio tests for the association between local rare variants and gene expression changes of three types: outlier genes expression, differential expression between parent and child, and alternative isoforms present in the child but not the parent.

### Libraries used

Python 3.7.3:

* pandas 1.0.0
* statsmodels 0.11.0

### outlier_log_odds.py

Calculates log odds ratios for the presence of rare variants near genes with outlier expression.

### diff_exp_log_odds.py

Calculates log odds ratios for the presence of rare variants near genes with differential expression between parent-child pairs.

### alt_splicing_log_odds.py

Calculates log odds ratios for the presence of rare variants near genes with unique isoforms between parent-child pairs.
