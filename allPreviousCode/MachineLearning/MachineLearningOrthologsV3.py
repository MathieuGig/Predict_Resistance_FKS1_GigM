# importing modules and packages
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV

import glob
import os

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)
########################################################################################################################
########################################################################################################################

# importing data & Concatenate all dataframes
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
                ])[['nt_seq']].nunique().reset_index()

shortMaster = master[['strain', 'locus', 'pool_type', 'seq_type', 'compound', 'median_s', 'alt_aa', 'aa_pos']]
shortMaster = shortMaster.loc[shortMaster['compound'].isin(['anidulafungin', 'caspofungin', 'micafungin'])]

########################################################################################################################

# Amino acid properties
AAproperties = pd.read_table('all_indices_final_table_propensity.txt')
AAproperties.rename(columns={'Aminoacid.1.letter': 'alt_aa'}, inplace=True)

########################################################################################################################

# Make wild-type dataframe.
WT_master = shortMaster.loc[(shortMaster['seq_type'] == 'WT') & (shortMaster['pool_type'] == 'single')]

# Merge with aa properties.
WT_merged = pd.merge(left=WT_master, right=AAproperties, how='inner', indicator='location', suffixes=(None, '_singles'),
                  on='alt_aa')
WT_merged = WT_merged.drop(columns=['seq_type', 'pool_type', 'location'])
WT_merged = WT_merged.sort_values(by=['strain', 'locus', 'compound', 'aa_pos'])

WT_sequences = WT_merged.groupby(['strain', 'locus', 'compound', 'median_s']).agg(''.join)
WT_properties = WT_merged.groupby(['strain', 'locus', 'compound', 'median_s']).median()
WT = pd.concat([WT_sequences, WT_properties], axis=1).reset_index()

########################################################################################################################

# Make orthologs dataframe.
Ortho_master = shortMaster.loc[shortMaster['seq_type'] == 'ortho']

# Merge with aa properties
Ortho_merged = pd.merge(left=Ortho_master, right=AAproperties, how='inner', indicator='location', suffixes=(None, '_singles'),
                  on='alt_aa')
Ortho_merged = Ortho_merged.drop(columns=['seq_type', 'pool_type', 'location'])
Ortho_merged = Ortho_merged.sort_values(by=['strain', 'locus', 'compound', 'aa_pos'])

Ortho_sequences = Ortho_merged.groupby(['strain', 'locus', 'compound', 'median_s']).agg(''.join)
Ortho_properties = Ortho_merged.groupby(['strain', 'locus', 'compound', 'median_s']).median()
Ortho = pd.concat([Ortho_sequences, Ortho_properties], axis=1).reset_index()

########################################################################################################################

# Make delta dataframe.