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

single = master[master.seq_type == 'single'][['strain','locus','compound','aa_pos','alt_codons','median_s']]
single.rename(columns={'median_s':'single_s'}, inplace=True)

print(single)

mergedf = pd.merge(left=master[master.seq_type.isin(['ortho'])],
                   right=single, how='inner',
                   on=['strain','locus','compound','aa_pos','alt_codons'])

merge_long = mergedf.pivot_table(index = ['strain', 'locus', 'compound','nt_seq','seq_type','median_s'],
                                 columns = 'mutated_codon', values = ['aa_pos', 'alt_codons','single_s'],
                                 aggfunc = 'first').reset_index()
merge_long = merge_long.fillna(0)

merge_long['additivity'] = merge_long['single_s', 1.0] + merge_long['single_s', 2.0] + merge_long['single_s', 3.0] + \
                           merge_long['single_s', 4.0] + merge_long['single_s', 5.0] + merge_long['single_s', 6.0] + \
                           merge_long['single_s', 7.0] + merge_long['single_s', 8.0] + merge_long['single_s', 9.0]

print(merge_long)



# Graphs

merge_long = merge_long.loc[merge_long['strain'] == 'BY4741']

compounds = ['caspofungin', 'micafungin', 'anidulafungin', 'none']
SortLocus = ['FKS1-HS1', 'FKS1-HS2', 'FKS2-HS1', 'FKS2-HS2']

graphdf = merge_long[(merge_long.compound.isin(compounds))]

fig = sns.relplot(graphdf, x='additivity', y='median_s',
                  row='locus', col='compound', col_order=compounds, row_order=SortLocus,
                  height=2.5, aspect=1)

fig.set(xlim=(-5, 10), ylim=(-5, 5))

xrange = np.linspace(-5, 10, 16)

for l in range(len(SortLocus)):
    for c in range(len(compounds)):
        fig.axes[l][c].axline((-5, -5), slope=1, ls='--', c='grey', zorder=-100)
        fig.axes[l][c].fill_between(xrange, xrange - 2, xrange + 2, color='grey', alpha=.4)

fig.set_titles(row_template='{row_name}', col_template='{col_name}')
fig.set(xlabel='Additivity', ylabel='Ortholog Sel_coef')

fig.fig.subplots_adjust(top=.88)
fig.fig.suptitle('BY4741')

plt.show()