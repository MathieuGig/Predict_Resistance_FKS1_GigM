# importing modules and packages
from random import randint

import pandas as pd
import numpy as np
np.bool = np.bool_
np.int = np.int_
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, accuracy_score, confusion_matrix, \
    ConfusionMatrixDisplay, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from scipy.stats import randint

import glob
import os

import shap # v0.39.0
shap.initjs()

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

########################################################################################################################

# VIP. Very important parameter
drug = 'caspofungin'

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

shortMaster = master[['strain', 'locus', 'pool_type', 'seq_type', 'compound', 'median_s', 'aa_seq', 'Nham_aa', 'alt_aa', 'aa_pos']]
shortMaster = shortMaster.loc[shortMaster['compound'].isin(['anidulafungin', 'caspofungin', 'micafungin'])]

########################################################################################################################

# Amino acid properties
AAproperties = pd.read_table('all_indices_final_table_propensity.txt')
AAproperties.rename(columns={'Aminoacid.1.letter': 'alt_aa'}, inplace=True)

########################################################################################################################

# Single mutations dataframe
Single_master = shortMaster.loc[(shortMaster['seq_type'] == 'single') & (shortMaster['pool_type'] == 'single') & (shortMaster['strain'] == 'BY4741') & (shortMaster['locus'] == 'FKS1-HS1') & (shortMaster['compound'] == drug)]
Single_master['Resistance'] = np.where(Single_master['median_s'] >= 1, 'resistant', 'susceptible')

Single_merged = pd.merge(left=Single_master, right=AAproperties, how='inner', indicator='location', suffixes=(None, '_singles'),
                  on='alt_aa')

X_train = Single_merged.drop(columns=['strain', 'locus', 'seq_type', 'aa_seq', 'pool_type', 'location', 'compound', 'median_s', 'Nham_aa', 'alt_aa', 'Resistance'])
y_train = Single_merged['Resistance']

########################################################################################################################

# Density plot of Singles
#sns.histplot(data=Single_master, x='median_s', element='step')
#plt.show()

########################################################################################################################

# Orthologs dataframe
Ortho_master = shortMaster.loc[(shortMaster['seq_type'] == 'ortho') & (shortMaster['pool_type'] == 'single') & (shortMaster['strain'] == 'BY4741') & (shortMaster['locus'] == 'FKS1-HS1') & (shortMaster['compound'] == drug)]
Ortho_master['Resistance'] = np.where(Ortho_master['median_s'] >= 1, 'resistant', 'susceptible')

Ortho_merged = pd.merge(left=Ortho_master, right=AAproperties, how='inner', indicator='location', suffixes=(None, '_singles'),
                  on='alt_aa')

X_test = Ortho_merged.drop(columns=['strain', 'locus', 'seq_type', 'aa_seq', 'pool_type', 'location', 'compound', 'median_s', 'Nham_aa', 'alt_aa', 'Resistance'])
y_test = Ortho_merged['Resistance']

########################################################################################################################

# Density plot of Orthologs
#sns.histplot(data=Ortho_master, x='median_s', element='step')
#plt.show()

#sns.histplot(data=Ortho_master, x='aa_pos', element='step')
#plt.show()

########################################################################################################################

#rf = RandomForestClassifier()
#rf.fit(X_train, y_train)
#y_pred = rf.predict(X_test)
#accuracy = accuracy_score(y_test, y_pred)
#print("Accuracy:", accuracy)
# Create the confusion matrix
#cm = confusion_matrix(y_test, y_pred)
#ConfusionMatrixDisplay(confusion_matrix=cm).plot();
#plt.show()


param_dist = {'n_estimators': randint(50,500),
              'max_depth': randint(4,40)}

# Create a random forest classifier
rf = RandomForestClassifier()

# Use random search to find the best hyperparameters
rand_search = RandomizedSearchCV(rf,
                                 param_distributions = param_dist,
                                 n_iter=20,
                                 cv=5, scoring='roc_auc')

# Fit the random search object to the data
rand_search.fit(X_train, y_train)

# Create a variable for the best model
best_rf = rand_search.best_estimator_

# Print the best hyperparameters
print('Best hyperparameters:',  rand_search.best_params_)

# Generate predictions with the best model
y_pred = best_rf.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)

# Create the confusion matrix
cm = confusion_matrix(y_test, y_pred)
ConfusionMatrixDisplay(confusion_matrix=cm).plot();
plt.show()
