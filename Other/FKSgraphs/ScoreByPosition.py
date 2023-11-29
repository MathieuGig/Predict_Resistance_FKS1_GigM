# importing modules and packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import glob
import os

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

# importing data. Concatenate all dataframes
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

    if pool_type != 'NaN':
        list_df.append(df)

master = pd.concat(list_df, ignore_index=True)
master.groupby(['strain', 'locus', 'compound'  # Per input
                ])[['aa_seq']].nunique().reset_index()

shortMaster = master[['strain', 'locus', 'pool_type', 'compound', 'median_s', 'alt_aa', 'aa_pos', 'WT', 'aa_seq']]

#strain = 'BY4741'
strain = 'R1158'

Condition = shortMaster.loc[shortMaster['strain'] == strain]

#Data = Condition.groupby(['strain', 'locus', 'compound'])[['aa_pos', 'median_s']].agg(aa_pos=('aa_pos', 'first'), median_s=('median_s', 'median')).reset_index()
Data = Condition.groupby(['strain', 'locus', 'compound', 'aa_pos'])['median_s'].agg('median').reset_index()

compounds = ['caspofungin', 'micafungin', 'anidulafungin', 'none']
#SortLocus = ['FKS1-HS1', 'FKS1-HS2']
SortLocus = ['FKS1-HS1', 'FKS1-HS2', 'FKS2-HS1', 'FKS2-HS2']


#### Make Graph

graphdf = Data[(Data.compound.isin(compounds))]

fig = sns.relplot(graphdf, x='aa_pos', y='median_s',
                  row='locus', col='compound', col_order=compounds, row_order=SortLocus,
                  height=2.5, aspect=1)

fig.set_titles(row_template='{row_name}', col_template='{col_name}')
fig.set(xlabel='mutated position', ylabel='Median Sel_coef')

fig.fig.subplots_adjust(top=.88)
fig.fig.suptitle(strain)

plt.show()