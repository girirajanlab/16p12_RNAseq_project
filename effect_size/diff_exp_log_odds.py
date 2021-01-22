import pandas as pd
from statsmodels.stats.contingency_tables import Table2x2

# load in outlier and rare variant data

df = pd.read_csv('../data/gene_level_expression.tsv', sep='\t')

# define rare variant classes

variants = list(df.columns)[7:]

# Create dataframe to populate

stats = pd.DataFrame(index=variants)

stats['rare_variant_and_exp_change'] = 0.
stats['rare_variant_and_not_exp_change'] = 0.
stats['not_rare_variant_and_exp_change'] = 0.
stats['not_rare_variant_and_not_exp_change'] = 0.
stats['odds_ratio'] = 0.
stats['log_odds_ratio'] = 0.
stats['log_odds_ratio_p_value'] = 0.
stats['log_odds_ratio_conf_lower'] = 0.
stats['log_odds_ratio_conf_upper'] = 0.

# Perform odds ratio test

for variant in variants:
    rare_var_and_exp_change        = df[(df[variant] != 0) & (~(df['Diff. expression'].isna()))].shape[0]
    rare_var_and_not_exp_change    = df[(df[variant] != 0) & (df['Diff. expression'].isna())].shape[0]
    no_rare_var_and_exp_change     = df[(df[variant] == 0) & (~(df['Diff. expression'].isna()))].shape[0]
    no_rare_var_and_not_exp_change = df[(df[variant] == 0) & (df['Diff. expression'].isna())].shape[0]

    stats.at[variant, 'rare_variant_and_exp_change']        = rare_var_and_exp_change
    stats.at[variant, 'rare_variant_and_not_exp_change']    = rare_var_and_not_exp_change 
    stats.at[variant, 'not_rare_variant_and_exp_change']     = no_rare_var_and_exp_change
    stats.at[variant, 'not_rare_variant_and_not_exp_change'] = no_rare_var_and_not_exp_change

    cont = Table2x2([[rare_var_and_exp_change, rare_var_and_not_exp_change], [no_rare_var_and_exp_change, no_rare_var_and_not_exp_change]])


    stats.at[variant, 'odds_ratio'] = cont.oddsratio
    stats.at[variant, 'log_odds_ratio'] = cont.log_oddsratio_confint()[0]
    stats.at[variant, 'log_odds_ratio_p_value'] = cont.log_oddsratio
    stats.at[variant, 'log_odds_ratio_conf_lower'] = cont.log_oddsratio_confint()[1]
    stats.at[variant, 'log_odds_ratio_conf_upper'] = cont.oddsratio_pvalue()


# save table
stats.index.name = 'Variant'
stats.to_csv('diff_exp_log_odds.tsv', sep='\t')