# importing modules and packages
import pandas as pd
import numpy as np
np.bool = np.bool_
np.int = np.int_
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

    X = merged.drop(columns=['strain', 'locus', 'compound', 'median_s', 'location'])
    X = X.drop(columns=['alt_aa'])
    y = merged['median_s']

    # creating train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=101)

    # Gridsearch & Random Forest
    grid = {
        'n_estimators': [200, 300, 400, 500, 600],
        'max_features': ['sqrt', 'log2', None],
        'max_depth': [9, 10, 11, 12, 13, None],
        'random_state': [18]
    }

    CV_rf = GridSearchCV(estimator=RandomForestRegressor(), param_grid=grid, n_jobs=-1)
    CV_rf.fit(X_train, y_train)
    CV_rf_pred = CV_rf.predict(X_test)
    print(CV_rf.best_params_)
    print('RandomForest r2 score : ', r2_score(y_test, CV_rf_pred))
    print('RandomForest MSE : ', mean_squared_error(y_test, CV_rf_pred))
    print('RandomForest MAE : ', mean_absolute_error(y_test, CV_rf_pred))

    plt.figure(figsize=(25, 12))
    feat_importances = pd.Series(CV_rf.best_estimator_.feature_importances_, index=X_train.columns)
    feat_importances.nlargest(5).plot(kind='barh')
    plt.title(
        f'Random Forest Regression feature importance : {strain} {locus} {drug} \n r2 score : {r2_score(y_test, CV_rf_pred):.3f}')
    plt.xlabel('coefficients')
    plt.ylabel('feature')
    plt.savefig(f'ML_pos1feature_{strain}_{locus}_{drug}.png')


GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'caspofungin')
GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'micafungin')
GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'anidulafungin')