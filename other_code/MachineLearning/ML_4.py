# importing modules and packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler

import glob
import os

#import shap # v0.39.0
#shap.initjs()

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
ConditionMaster = shortMaster.loc[shortMaster['strain'] == 'BY4741']
ConditionMaster = ConditionMaster.loc[ConditionMaster['locus'] == 'FKS1-HS1']
ConditionMaster = ConditionMaster.loc[ConditionMaster['pool_type'] == 'single']
ConditionMaster = ConditionMaster.loc[ConditionMaster['compound'] == 'caspofungin']
ConditionMaster = ConditionMaster.loc[ConditionMaster['Nmut_codons'] == 1.0]

ModelData = ConditionMaster.groupby(['strain', 'locus', 'compound', 'aa_seq', 'alt_aa'])[['median_s', 'aa_pos']].median().reset_index()

# Amino acid properties
AAproperties = pd.read_table('all_indices_final_table_propensity.txt')
AAproperties.rename(columns={'Aminoacid.1.letter':'alt_aa'}, inplace=True)

# Merge data frames.
merged = pd.merge(left=ModelData, right=AAproperties, how='inner', indicator='location', suffixes=(None, '_singles'), on='alt_aa')
print(merged)

#Creating feature variables
X = merged.drop(columns=['strain', 'locus', 'compound', 'median_s', 'alt_aa', 'aa_seq', 'location'])
y = merged['median_s']

# creating train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=101)

# Standardization
scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_test_s = scaler.transform(X_test)

# creating a regression model
model = LinearRegression()

# fitting the model
model.fit(X_train_s, y_train)


# shap
#explainer = shap.Explainer(model)
#shap_test = explainer(X_test)
#shap_df = pd.DataFrame(shap_test.values, columns=shap_test.feature_names, index=X_test.index)

#columns = shap_df.apply(np.abs).mean()\
#                 .sort_values(ascending=False).indexfig, ax = plt.subplots(1, 2, figsize=(11,4))
#sns.barplot(data=shap_df[columns].apply(np.abs), orient='h',
#            ax=ax[0])
#ax[0].set_title("Mean absolute shap value")
#sns.boxplot(data=shap_df[columns], orient='h', ax=ax[1])
#ax[1].set_title("Distribution of shap values");

# making predictions
predictions = model.predict(X_test_s)

# model evaluation
score = r2_score(y_test, predictions)
print('r2 score : ', score)

ms_error = mean_squared_error(y_test, predictions)
print('mean_squared_error : ', ms_error)

ma_error = mean_absolute_error(y_test, predictions)
print('mean_absolute_error : ', ma_error)

# graph
plt.figure(figsize=(10,7))
feat_importances = pd.Series(model.coef_, index= X_train.columns)
feat_importances.nlargest(65).plot(kind='barh')
plt.title('by4741 fks1 hs1 single caspo')
plt.xlabel('feature importance')
plt.ylabel('feature')
plt.show()


train_dataset = X_train.copy()
train_dataset.insert(0, "median_s", y_train)
_ = sns.pairplot(train_dataset[['retention_coefficient_tfa_browne', 'retention_coefficient_hfba_browne', 'median_s', 'aa_pos']],
                 kind="reg",
                 diag_kind="kde",
                 plot_kws={"scatter_kws": {"alpha": 0.1}},)
plt.show()