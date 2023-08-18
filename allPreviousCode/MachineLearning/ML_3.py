# importing modules and packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn import preprocessing

from sklearn.metrics import RocCurveDisplay

import glob
import os

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

# importing data

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

    if pool_type != 'NaN':
        list_df.append(df)

master = pd.concat(list_df, ignore_index=True)
master.groupby(['strain', 'locus', 'compound'  # Per input
                ])[['nt_seq']].nunique().reset_index()

shortMaster = master[['strain', 'locus', 'pool_type', 'compound', 'median_s', 'alt_aa', 'aa_pos', 'Nham_aa', 'Nmut_codons', 'WT', 'aa_seq']]

# Choose strain, locus, compound, and pool type.
NewMaster = shortMaster.loc[shortMaster['strain'] == 'BY4741']
NewMaster = NewMaster.loc[NewMaster['locus'] == 'FKS1-HS1']
NewMaster = NewMaster.loc[NewMaster['pool_type'] == 'single']
NewMaster = NewMaster.loc[NewMaster['compound'] == 'caspofungin']


#print(NewMaster)
#iam = NewMaster[['aa_seq', 'median_s']]
#miam = miam.groupby('aa_seq').median()
#miam.rename(columns={'median_s':'median_median_s'}, inplace=True)
#medianMerge = pd.merge(left=NewMaster, right=miam, how='inner', indicator='ok', suffixes=(None, 'ok'), on='aa_seq')
#print(medianMerge)
#NewMedianMerge = medianMerge.drop(columns=['median_s', 'ok'])
#WOWMedianMerge = NewMedianMerge.drop_duplicates()
#print(WOWMedianMerge)

NewNewMaster = NewMaster.groupby(['strain', 'locus', 'compound', 'aa_seq', 'alt_aa'])[['median_s']].median().reset_index()
print(NewNewMaster)

# Amino acid properties
AAproperties = pd.read_table('all_indices_final_table_propensity.txt')
AAproperties.rename(columns={'Aminoacid.1.letter':'alt_aa'}, inplace=True)

# Merge data frames.
merged = pd.merge(left=NewNewMaster, right=AAproperties, how='inner', indicator='location', suffixes=(None, '_singles'), on='alt_aa')

#Creating feature variables
#X = merged.drop(columns=['strain', 'locus', 'pool_type', 'compound', 'median_s', 'alt_aa', 'WT', 'Nham_aa', 'Nmut_codons', 'aa_seq', 'location'])
X = merged.drop(columns=['strain', 'locus', 'compound', 'median_s', 'alt_aa', 'aa_seq', 'location'])
y = merged['median_s']

# creating train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=101)

# creating a regression model
model = LinearRegression()

# fitting the model
model.fit(X_train, y_train)

# making predictions
predictions = model.predict(X_test)

# model evaluation
print('mean_squared_error : ', mean_squared_error(y_test, predictions))
print('mean_absolute_error : ', mean_absolute_error(y_test, predictions))

# graph
plt.figure(figsize=(10,7))
feat_importances = pd.Series(model.coef_, index= X_train.columns)
feat_importances.nlargest(7).plot(kind='barh')
plt.title('by4741 fks1 hs1 single caspo')
plt.xlabel('feature importance')
plt.ylabel('feature')
plt.show()