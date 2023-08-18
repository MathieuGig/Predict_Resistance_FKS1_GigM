# importing modules and packages
import pandas as pd
import numpy as np
np.bool = np.bool_
np.int = np.int_
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV

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

shortMaster = master[['strain', 'locus', 'pool_type', 'compound', 'median_s', 'alt_aa', 'aa_pos', 'Nham_aa', 'Nmut_codons', 'WT']]

def GenerateFeatureGraphML(strain, locus, drug):

    ConditionMaster = shortMaster.loc[shortMaster['strain'] == strain]
    ConditionMaster = ConditionMaster.loc[ConditionMaster['locus'] == locus]
    ConditionMaster = ConditionMaster.loc[ConditionMaster['pool_type'] == 'single']
    ConditionMaster = ConditionMaster.loc[ConditionMaster['compound'] == drug]
    ConditionMaster = ConditionMaster.loc[ConditionMaster['Nmut_codons'] == 1.0]

    ModelData = ConditionMaster.groupby(['strain', 'locus', 'compound', 'aa_pos', 'alt_aa'])[['median_s']].median().reset_index()

    # Amino acid properties
    AAproperties = pd.read_table('all_indices_final_table_propensity.txt')
    AAproperties.rename(columns={'Aminoacid.1.letter':'alt_aa'}, inplace=True)

    # Merge data frames.
    merged = pd.merge(left=ModelData, right=AAproperties, how='inner', indicator='location', suffixes=(None, '_singles'), on='alt_aa')

    #Creating feature variables
    X = merged.drop(columns=['strain', 'locus', 'compound', 'median_s', 'location'])
    X["pos&alt"] = X['aa_pos'].astype(str) + "_" + X["alt_aa"]
    X = X.drop(columns=['aa_pos', 'alt_aa'])
    Mutations = pd.get_dummies(X['pos&alt'], prefix='mutation')
    X = pd.concat([X, Mutations], axis=1)
    X = X.drop(columns=['pos&alt'])
    print(X)
    """
    y = merged['median_s']

    # creating train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=101)

    # Gridsearch & Random Forest
    grid = {
        'n_estimators': [500, 600, 700, 800, 900, 1000],
        'max_features': [None],
        'max_depth': [2, 3, 4, 5, 6, 7, 8, 9, 10, None],
        'random_state': [18]
    }

    CV_rf = GridSearchCV(estimator=RandomForestRegressor(), param_grid=grid, n_jobs=-1)
    CV_rf.fit(X_train, y_train)
    CV_rf_pred = CV_rf.predict(X_test)
    print(CV_rf.best_params_)

    plt.figure(figsize=(25, 12))
    feat_importances = pd.Series(CV_rf.best_estimator_.feature_importances_, index=X_train.columns)
    feat_importances.nlargest(5).plot(kind='barh')
    plt.title(f'Random Forest Regression feature importance : {strain} {locus} {drug} \n r2 score : {r2_score(y_test, CV_rf_pred):.3f}')
    plt.xlabel('coefficients')
    plt.ylabel('feature')

    plt.savefig(f'ML_DMSDummy_{strain}_{locus}_{drug}.png')
    """


GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'caspofungin')