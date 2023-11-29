import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import glob
import os
import numpy as np

# Voir pandas.tuple


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

#On aime lui
#per_mut = master[master.seq_type.isin(['single', 'ortho'])].pivot_table(
    #index=['strain', 'locus', 'compound', 'nt_seq', 'seq_type', 'alt_codons'],
    #columns='aa_pos', values='median_s',
    #aggfunc='first').reset_index()

pd.set_option('mode.chained_assignment', None)


Orthologs = master[master.seq_type.isin(['ortho'])].pivot_table(
    index=['strain', 'locus', 'compound', 'nt_seq', 'seq_type', 'median_s', 'alt_codons'],
    columns='aa_pos', values='mutated_codon',
    aggfunc='first').reset_index()

Singles = master[master.seq_type.isin(['single'])].pivot_table(
    index=['strain', 'locus', 'compound', 'nt_seq', 'seq_type', 'alt_codons'],
    columns='aa_pos', values='median_s',
    aggfunc='first').reset_index()
#Singles.rename(columns={'median_s': 's_single'}, inplace=True)

Merged = pd.merge(left=Orthologs, right=Singles, how='inner', indicator='location', suffixes=(None, '_singles'),
                          on=['strain', 'locus', 'compound', 'alt_codons'])

print(Merged)