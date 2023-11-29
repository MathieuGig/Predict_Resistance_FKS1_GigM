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

np.bool = np.bool_
np.int = np.int_
import shap # v0.39.0
shap.initjs()

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

shortMaster = master[['strain', 'locus', 'pool_type', 'seq_type', 'compound', 'median_s', 'aa_seq', 'alt_aa', 'aa_pos']]
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
WT_merged = WT_merged.drop(columns=['seq_type', 'aa_seq', 'pool_type', 'location'])
WT_merged = WT_merged.sort_values(by=['strain', 'locus', 'compound', 'aa_pos'])

WT_sequences = WT_merged.groupby(['strain', 'locus', 'compound', 'median_s']).agg(''.join)
WT_properties = WT_merged.groupby(['strain', 'locus', 'compound', 'median_s']).median()
WT = pd.concat([WT_sequences, WT_properties], axis=1).reset_index()

########################################################################################################################

# Make orthologs dataframe.
Ortho_master = shortMaster.loc[shortMaster['seq_type'] == 'ortho']
orthoSeq_list = list(Ortho_master['aa_seq'])
output_list = list(map(list, orthoSeq_list))
Ortho_master['aa_seq'] = output_list
Ortho_master = Ortho_master.explode('aa_seq')

# Merge with aa properties
AAproperties.rename(columns={'alt_aa': 'aa_seq'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_master, right=AAproperties, how='inner', indicator='location', suffixes=(None, '_merged'),
                  on='aa_seq')
Ortho_merged = Ortho_merged.drop(columns=['seq_type', 'pool_type', 'location'])
Ortho_merged = Ortho_merged.sort_values(by=['strain', 'locus', 'compound', 'alt_aa', 'aa_pos'])

Ortho_sequences = Ortho_merged.groupby(['strain', 'locus', 'compound', 'median_s', 'aa_pos', 'alt_aa']).agg(''.join)
Ortho_properties = Ortho_merged.groupby(['strain', 'locus', 'compound', 'median_s', 'aa_pos', 'alt_aa']).median()
Ortho = pd.concat([Ortho_sequences, Ortho_properties], axis=1).reset_index()
Ortho = Ortho.drop(columns=['aa_pos', 'alt_aa'])
Ortho = Ortho.groupby(['strain', 'locus', 'compound', 'median_s', 'aa_seq']).median().reset_index()

########################################################################################################################

# Make delta dataframe.
WT = WT.drop(columns=['alt_aa', 'aa_pos'])
Ortho = Ortho.drop(columns='aa_seq')

def GenerateFeatureGraphML(strain, locus, drug):
    print(f'{strain}_{locus}_{drug}')

    Condition_WT = WT.loc[(WT['strain'] == strain) & (WT['locus'] == locus) & (WT['compound'] == drug)]
    Condition_WT = Condition_WT.drop(columns=['median_s'])

    Condition_Ortho = Ortho.loc[(Ortho['strain'] == strain) & (Ortho['locus'] == locus) & (Ortho['compound'] == drug)]
    y = Condition_Ortho['median_s']
    Condition_Ortho = Condition_Ortho.drop(columns=['median_s'])

    Delta = Condition_Ortho.merge(Condition_WT, how='left', on=['strain', 'locus', 'compound'], suffixes=('', '_WT'))
    Delta = Delta.drop(columns=['strain', 'locus', 'compound'])
    for i in Delta.columns:
        if '_WT' not in i:
            Delta[f'delta_{i}'] = Delta[f'{i}'] - Delta[f'{i}_WT']
            Delta = Delta.drop(columns=[f'{i}', f'{i}_WT'])

    #print(f'Number of rows : {len(Delta)}')
    ########################################################################################################################

    # Machine Learning

    # creating train and test sets
    X_train, X_test, y_train, y_test = train_test_split(Delta, y, test_size=0.3, random_state=101)

    # Gridsearch & Random Forest
    grid = {'n_estimators': [600, 700, 800, 900, 1000],
            'max_features': ['sqrt', 'log2', None],
            'max_depth': [5, 6, 7, 8, 9, 10, None],
            'random_state': [18]
    }

    CV_rf = GridSearchCV(estimator=RandomForestRegressor(), param_grid=grid, n_jobs=-1)
    CV_rf.fit(X_train, y_train)
    CV_rf_pred = CV_rf.predict(X_test)
    print(CV_rf.best_params_)
    print('RandomForest r2 score : ', r2_score(y_test, CV_rf_pred))
    print('RandomForest MSE : ', mean_squared_error(y_test, CV_rf_pred))
    print('RandomForest MAE : ', mean_absolute_error(y_test, CV_rf_pred))


    # Standardization
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    # create a dummy regressor
    dummy_reg = DummyRegressor(strategy='mean')
    # fit it on the training set
    dummy_reg.fit(X_train_s, y_train)
    # make predictions on the test set
    dummy_pred = dummy_reg.predict(X_test_s)
    print('Dummy r2 score : ', r2_score(y_test, dummy_pred))
    #print('Dummy MSE : ', mean_squared_error(y_test, dummy_pred))
    #print('Dummy MAE : ', mean_absolute_error(y_test, dummy_pred))


    #plt.figure(figsize=(25, 25))
    #feat_importances = pd.Series(CV_rf.best_estimator_.feature_importances_, index=X_train.columns)
    #feat_importances.nlargest(5).plot(kind='barh')
    plt.title(
        f'Random Forest Regression feature importance : {strain} {locus} {drug} \n model r2 score : {r2_score(y_test, CV_rf_pred):.3f}',
    fontsize=20)
    #plt.xlabel('coefficients')
    #plt.ylabel('feature')

    # shap
    model = CV_rf.best_estimator_
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)
    shap.summary_plot(shap_values, X_test, show=False, plot_size=(16,8))
    plt.savefig(f'ML_shap_{strain}_{locus}_{drug}.png')

GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'caspofungin')
#GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'micafungin')
#GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'anidulafungin')

#GenerateFeatureGraphML('R1158', 'FKS1-HS1', 'caspofungin')
#GenerateFeatureGraphML('R1158', 'FKS1-HS1', 'micafungin')
#GenerateFeatureGraphML('R1158', 'FKS1-HS1', 'anidulafungin')
#GenerateFeatureGraphML('R1158', 'FKS1-HS2', 'caspofungin')
#GenerateFeatureGraphML('R1158', 'FKS1-HS2', 'micafungin')
#GenerateFeatureGraphML('R1158', 'FKS1-HS2', 'anidulafungin')

#GenerateFeatureGraphML('R1158', 'FKS2-HS1', 'caspofungin')
#GenerateFeatureGraphML('R1158', 'FKS2-HS1', 'micafungin')
#GenerateFeatureGraphML('R1158', 'FKS2-HS1', 'anidulafungin')
#GenerateFeatureGraphML('R1158', 'FKS2-HS2', 'caspofungin')
#GenerateFeatureGraphML('R1158', 'FKS2-HS2', 'micafungin')