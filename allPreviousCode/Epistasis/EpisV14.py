import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import glob
import os
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

# Concatenate all dataframes
list_df = []
for d in glob.glob('DMS-main/processed_data/*'):
    cset = os.path.basename(d).split('_')
    cset_name = '_'.join(cset)
    strain, locus, pool_type, compound = [cset[i] for i in (0, 1, 2, -1)]

    df = pd.read_csv(f'{d}/selcoeff_all_libraries.csv', index_col=0)
    df['strain'] = strain
    df['locus'] = locus
    df['pool_type'] = pool_type
    df['compound'] = compound

    if pool_type != 'double':
        list_df.append(df)

master = pd.concat(list_df, ignore_index=True)
master.groupby(['strain', 'locus', 'compound'  # Per input
                ])[['nt_seq']].nunique().reset_index()

pd.set_option('mode.chained_assignment', None)

Orthologs = master[master.seq_type.isin(['ortho'])].pivot_table(
    index=['strain', 'locus', 'compound', 'nt_seq', 'seq_type', 'median_s'],
    columns='aa_pos', values='alt_codons',
    aggfunc='first').reset_index()

Singles = master[master.seq_type.isin(['single'])].pivot_table(
    index=['strain', 'locus', 'compound', 'nt_seq', 'seq_type', 'median_s'],
    columns='aa_pos', values='alt_codons',
    aggfunc='first').reset_index()

Merged = pd.merge(left=Orthologs, right=Singles, how='inner', indicator='location', suffixes=(None, '_singles'),
                          on=['strain', 'locus', 'compound', 0.0])


#Try the loop!

#for i in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
    #per_mut_concat_merged = pd.merge(left=per_mut_merged_both, right=per_mut_single, how='outer', indicator=f'location_{i}', suffixes=(None, f'_{i}'),
                              #on=['strain', 'locus', 'compound', f'aa_pos_{i}', i])

    #per_mut_merged_both = per_mut_concat_merged.loc[per_mut_concat_merged[f'location_{i}'] == 'both']


print(Merged)