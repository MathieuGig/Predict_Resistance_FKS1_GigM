# importing modules and packages
import pandas as pd
import numpy as np
np.bool = np.bool_
np.int = np.int_
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, matthews_corrcoef, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV

import glob
import os

import shap # v0.39.0
shap.initjs()

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

########################################################################################################################

# VIP. Very important parameter. The antifungal drugs.
drug2 = 'caspofungin'

drug1 = 'anidulafungin'

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
AAproperties.rename(columns={'Aminoacid.1.letter': 'aa1'}, inplace=True)

########################################################################################################################

# Single mutations dataframe on drug 1
Drug1_master = shortMaster.loc[(shortMaster['seq_type'] == 'single') & (shortMaster['pool_type'] == 'single') & (shortMaster['strain'] == 'BY4741') & (shortMaster['locus'] == 'FKS1-HS1') & (shortMaster['compound'] == drug1)]
Drug1_master['Resistance'] = np.where(Drug1_master['median_s'] >= 1, 'resistant', 'susceptible')

Drug1_master[['aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9']] = Drug1_master['aa_seq'].apply(lambda x: pd.Series(list(x)))

Drug1_merged = pd.merge(left=Drug1_master, right=AAproperties, how='inner', indicator='location1', suffixes=(None, '_aa1'),
                  on='aa1')

AAproperties.rename(columns={'aa1': 'aa2'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator='location2', suffixes=(None, '_aa2'),
                  on='aa2')

AAproperties.rename(columns={'aa2': 'aa3'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator='location3', suffixes=(None, '_aa3'),
                  on='aa3')

AAproperties.rename(columns={'aa3': 'aa4'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator='location4', suffixes=(None, '_aa4'),
                  on='aa4')

AAproperties.rename(columns={'aa4': 'aa5'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator='location5', suffixes=(None, '_aa5'),
                  on='aa5')

AAproperties.rename(columns={'aa5': 'aa6'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator='location6', suffixes=(None, '_aa6'),
                  on='aa6')

AAproperties.rename(columns={'aa6': 'aa7'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator='location7', suffixes=(None, '_aa7'),
                  on='aa7')

AAproperties.rename(columns={'aa7': 'aa8'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator='location8', suffixes=(None, '_aa8'),
                  on='aa8')

AAproperties.rename(columns={'aa8': 'aa9'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator='location9', suffixes=(None, '_aa9'),
                  on='aa9')

X_train = Drug1_merged.drop(columns=['strain', 'locus', 'seq_type', 'aa_seq', 'pool_type', 'compound', 'median_s', 'Nham_aa', 'alt_aa', 'Resistance',
                                      'location1', 'location2', 'location3', 'location4', 'location5', 'location6', 'location7', 'location8', 'location9',
                                      'aa_pos', 'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9'])
y_train = Drug1_merged['Resistance']

########################################################################################################################

# Single mutations dataframe on drug 2
Drug2_master = shortMaster.loc[(shortMaster['seq_type'] == 'single') & (shortMaster['pool_type'] == 'single') & (shortMaster['strain'] == 'BY4741') & (shortMaster['locus'] == 'FKS1-HS1') & (shortMaster['compound'] == drug2)]
Drug2_master['Resistance'] = np.where(Drug2_master['median_s'] >= 1, 'resistant', 'susceptible')
Drug2_master = Drug2_master.drop(columns=['aa_pos', 'alt_aa'])
Drug2_master = Drug2_master.drop_duplicates()
Drug2_master = Drug2_master.reset_index(drop=True)

Drug2_master[['aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9']] = Drug2_master['aa_seq'].apply(lambda x: pd.Series(list(x)))

AAproperties.rename(columns={'aa9': 'aa1'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_master, right=AAproperties, how='inner', indicator='location1', suffixes=(None, '_aa1'),
                  on='aa1')

AAproperties.rename(columns={'aa1': 'aa2'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator='location2', suffixes=(None, '_aa2'),
                  on='aa2')

AAproperties.rename(columns={'aa2': 'aa3'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator='location3', suffixes=(None, '_aa3'),
                  on='aa3')

AAproperties.rename(columns={'aa3': 'aa4'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator='location4', suffixes=(None, '_aa4'),
                  on='aa4')

AAproperties.rename(columns={'aa4': 'aa5'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator='location5', suffixes=(None, '_aa5'),
                  on='aa5')

AAproperties.rename(columns={'aa5': 'aa6'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator='location6', suffixes=(None, '_aa6'),
                  on='aa6')

AAproperties.rename(columns={'aa6': 'aa7'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator='location7', suffixes=(None, '_aa7'),
                  on='aa7')

AAproperties.rename(columns={'aa7': 'aa8'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator='location8', suffixes=(None, '_aa8'),
                  on='aa8')

AAproperties.rename(columns={'aa8': 'aa9'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator='location9', suffixes=(None, '_aa9'),
                  on='aa9')

X_test = Drug2_merged.drop(columns=['strain', 'locus', 'seq_type', 'aa_seq', 'pool_type', 'compound', 'median_s', 'Nham_aa', 'Resistance',
                                      'location1', 'location2', 'location3', 'location4', 'location5', 'location6', 'location7', 'location8', 'location9',
                                    'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9'])
y_test = Drug2_merged['Resistance']

########################################################################################################################
## Machine Learning

# Gridsearch & Random Forest
#grid = {'n_estimators': [75, 100, 125, 150, 200, 250],
#        'max_features': ['sqrt', 'log2', None],
#        'max_depth': [5, 6, 7, None],
#        'random_state': [18]
#}

#CV_rf = GridSearchCV(estimator=RandomForestClassifier(), param_grid=grid, n_jobs=-1, cv=10)
#CV_rf.fit(X_train, y_train)

#rf = CV_rf.best_estimator_
#print(CV_rf.best_params_)

#y_pred = rf.predict(X_test)

rf = RandomForestClassifier(n_estimators= 250, max_features= 'sqrt', max_depth= 5, random_state= 18)

rf.fit(X_train, y_train)
y_pred = rf.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)

# Matthews correlation coefficient
mat = matthews_corrcoef(y_test, y_pred)
print("Mat:", mat)
#print(classification_report(y_test, y_pred))

# Create the confusion matrix
cmatrix = confusion_matrix(y_test, y_pred)

# Display the confusion matrix using seaborn heatmap
plt.figure(figsize=(10, 10))
sns.set(font_scale=3)
sns.heatmap(cmatrix, annot=True, fmt='d', cmap='Blues', xticklabels=['resistant', 'susceptible'], yticklabels=['resistant', 'susceptible'], annot_kws={"size": 60})
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title(f'Caspofungin Confusion Matrix \n Accuracy : {accuracy:.3f}')
plt.tight_layout()
#plt.savefig(f'Conf_matrix.svg')
plt.show()