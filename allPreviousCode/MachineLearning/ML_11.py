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
    Positions = pd.get_dummies(X['aa_pos'], prefix='aa_position')
    #Mutations = pd.get_dummies(X['alt_aa'], prefix='alt_aa')
    X = X.drop(columns=['aa_pos', 'alt_aa'])
    X = pd.concat([X, Positions], axis=1)
    #X = pd.concat([X, Mutations], axis=1)
    y = merged['median_s']

    # creating train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=101)

    # Standardization
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    X_test_s_df = pd.DataFrame(X_test_s, columns=X_test.columns)

    # creating a regression model
    model = LinearRegression()

    # fitting the model
    model.fit(X_train_s, y_train)

    # making predictions
    predictions = model.predict(X_test_s)


    # model evaluation
    score = r2_score(y_test, predictions)
    print('r2 score : ', score)

    ms_error = mean_squared_error(y_test, predictions)
    print('mean_squared_error : ', ms_error)

    ma_error = mean_absolute_error(y_test, predictions)
    print('mean_absolute_error : ', ma_error)

    # create a dummy regressor
    dummy_reg = DummyRegressor(strategy='mean')
    # fit it on the training set
    dummy_reg.fit(X_train_s, y_train)
    # make predictions on the test set
    dummy_pred = dummy_reg.predict(X_test_s)
    print('Dummy r2 score : ', r2_score(y_test, dummy_pred))
    print('Dummy MSE : ', mean_squared_error(y_test, dummy_pred))
    print('Dummy MAE : ', mean_absolute_error(y_test, dummy_pred))

    rf = RandomForestRegressor(n_estimators=900, max_features=None, max_depth=8, random_state=18)
    rf.fit(X_train_s, y_train)
    rf_pred = rf.predict(X_test_s)
    print('RandomForest r2 score : ', r2_score(y_test, rf_pred))
    print('RandomForest MSE : ', mean_squared_error(y_test, rf_pred))
    print('RandomForest MAE : ', mean_absolute_error(y_test, rf_pred))


    # Feature Importance Graph

    fig, axes = plt.subplots(1, 2, sharex=True, figsize=(25,10))
    fig.suptitle(f'{strain} {locus} {drug}')
    axes[0].set_title(f'Regression \n r2 score : {score}')
    axes[1].set_title(f'Random Forest Regression \n r2 score : {r2_score(y_test, rf_pred)}')


    feat_importances = pd.Series(model.coef_, index=X_train.columns)
    feat_importances.nlargest(5).plot(kind='barh', ax=axes[0])

    feat_importances2 = pd.Series(rf.feature_importances_, index=X_train.columns)
    feat_importances2.nlargest(5).plot(kind='barh', ax=axes[1])

    plt.savefig(f'ML_{strain}_{locus}_{drug}.png')

GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'caspofungin')
GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'micafungin')
GenerateFeatureGraphML('BY4741', 'FKS1-HS1', 'anidulafungin')

GenerateFeatureGraphML('BY4741', 'FKS1-HS2', 'caspofungin')
GenerateFeatureGraphML('BY4741', 'FKS1-HS2', 'micafungin')
GenerateFeatureGraphML('BY4741', 'FKS1-HS2', 'anidulafungin')

GenerateFeatureGraphML('R1158', 'FKS1-HS1', 'caspofungin')
GenerateFeatureGraphML('R1158', 'FKS1-HS1', 'micafungin')
GenerateFeatureGraphML('R1158', 'FKS1-HS1', 'anidulafungin')
GenerateFeatureGraphML('R1158', 'FKS1-HS2', 'caspofungin')
GenerateFeatureGraphML('R1158', 'FKS1-HS2', 'micafungin')
GenerateFeatureGraphML('R1158', 'FKS1-HS2', 'anidulafungin')

GenerateFeatureGraphML('R1158', 'FKS2-HS1', 'caspofungin')
GenerateFeatureGraphML('R1158', 'FKS2-HS1', 'micafungin')
GenerateFeatureGraphML('R1158', 'FKS2-HS1', 'anidulafungin')
GenerateFeatureGraphML('R1158', 'FKS2-HS2', 'caspofungin')
GenerateFeatureGraphML('R1158', 'FKS2-HS2', 'micafungin')
GenerateFeatureGraphML('R1158', 'FKS2-HS2', 'anidulafungin')