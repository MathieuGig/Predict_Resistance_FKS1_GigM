# importing modules and packages
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV

np.bool = np.bool_
np.int = np.int_
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, matthews_corrcoef, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import glob
import os

import shap # v0.39.0
shap.initjs()

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

########################################################################################################################

# VIP. Very important parameter. The Antifungal drug.
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
AAproperties.rename(columns={'Aminoacid.1.letter': 'aa1'}, inplace=True)

########################################################################################################################

## Standardization & PCA

# AAproperties[AAproperties.columns[1:]] = (AAproperties[AAproperties.columns[1:]]-AAproperties[AAproperties.columns[1:]].mean())/AAproperties[AAproperties.columns[1:]].std()
AAproperties[AAproperties.columns[1:]] = StandardScaler().fit_transform(AAproperties[AAproperties.columns[1:]])


pca = PCA(n_components=15)
wow = pca.fit_transform(AAproperties[AAproperties.columns[1:]])
print(wow)

########################################################################################################################
"""
# Single mutations dataframe
Single_master = shortMaster.loc[(shortMaster['seq_type'] == 'single') & (shortMaster['pool_type'] == 'single') & (shortMaster['strain'] == 'BY4741') & (shortMaster['locus'] == 'FKS1-HS1') & (shortMaster['compound'] == drug)]
Single_master['Resistance'] = np.where(Single_master['median_s'] >= 1, 'resistant', 'susceptible')

Single_master[['aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9']] = Single_master['aa_seq'].apply(lambda x: pd.Series(list(x)))

Single_merged = pd.merge(left=Single_master, right=AAproperties, how='inner', indicator='location1', suffixes=(None, '_aa1'),
                  on='aa1')

AAproperties.rename(columns={'aa1': 'aa2'}, inplace=True)
Single_merged = pd.merge(left=Single_merged, right=AAproperties, how='inner', indicator='location2', suffixes=(None, '_aa2'),
                  on='aa2')

AAproperties.rename(columns={'aa2': 'aa3'}, inplace=True)
Single_merged = pd.merge(left=Single_merged, right=AAproperties, how='inner', indicator='location3', suffixes=(None, '_aa3'),
                  on='aa3')

AAproperties.rename(columns={'aa3': 'aa4'}, inplace=True)
Single_merged = pd.merge(left=Single_merged, right=AAproperties, how='inner', indicator='location4', suffixes=(None, '_aa4'),
                  on='aa4')

AAproperties.rename(columns={'aa4': 'aa5'}, inplace=True)
Single_merged = pd.merge(left=Single_merged, right=AAproperties, how='inner', indicator='location5', suffixes=(None, '_aa5'),
                  on='aa5')

AAproperties.rename(columns={'aa5': 'aa6'}, inplace=True)
Single_merged = pd.merge(left=Single_merged, right=AAproperties, how='inner', indicator='location6', suffixes=(None, '_aa6'),
                  on='aa6')

AAproperties.rename(columns={'aa6': 'aa7'}, inplace=True)
Single_merged = pd.merge(left=Single_merged, right=AAproperties, how='inner', indicator='location7', suffixes=(None, '_aa7'),
                  on='aa7')

AAproperties.rename(columns={'aa7': 'aa8'}, inplace=True)
Single_merged = pd.merge(left=Single_merged, right=AAproperties, how='inner', indicator='location8', suffixes=(None, '_aa8'),
                  on='aa8')

AAproperties.rename(columns={'aa8': 'aa9'}, inplace=True)
Single_merged = pd.merge(left=Single_merged, right=AAproperties, how='inner', indicator='location9', suffixes=(None, '_aa9'),
                  on='aa9')

X_train = Single_merged.drop(columns=['strain', 'locus', 'seq_type', 'aa_seq', 'pool_type', 'compound', 'median_s', 'Nham_aa', 'alt_aa', 'Resistance',
                                      'location1', 'location2', 'location3', 'location4', 'location5', 'location6', 'location7', 'location8', 'location9',
                                      'aa_pos', 'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9'])
y_train = Single_merged['Resistance']

########################################################################################################################

# Density plot of Single mutations
sns.histplot(data=Single_master, x='median_s', element='step')
plt.title('Resistance Score Histogram', fontsize=24)
plt.xlabel('Resistance score', fontsize=20)
plt.ylabel('Count', fontsize=20)
plt.gca().tick_params(labelsize=18)
plt.tight_layout()
plt.savefig('SinglesDensity.svg')
plt.show()

########################################################################################################################

# Orthologs dataframe
Ortho_master = shortMaster.loc[(shortMaster['seq_type'] == 'ortho') & (shortMaster['pool_type'] == 'single') & (shortMaster['strain'] == 'BY4741') & (shortMaster['locus'] == 'FKS1-HS1') & (shortMaster['compound'] == drug)]
Ortho_master['Resistance'] = np.where(Ortho_master['median_s'] >= 1, 'resistant', 'susceptible')
Ortho_master = Ortho_master.drop(columns=['aa_pos', 'alt_aa'])
Ortho_master = Ortho_master.drop_duplicates()

Ortho_master[['aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9']] = Ortho_master['aa_seq'].apply(lambda x: pd.Series(list(x)))

AAproperties.rename(columns={'aa9': 'aa1'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_master, right=AAproperties, how='inner', indicator='location1', suffixes=(None, '_aa1'),
                  on='aa1')

AAproperties.rename(columns={'aa1': 'aa2'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_merged, right=AAproperties, how='inner', indicator='location2', suffixes=(None, '_aa2'),
                  on='aa2')

AAproperties.rename(columns={'aa2': 'aa3'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_merged, right=AAproperties, how='inner', indicator='location3', suffixes=(None, '_aa3'),
                  on='aa3')

AAproperties.rename(columns={'aa3': 'aa4'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_merged, right=AAproperties, how='inner', indicator='location4', suffixes=(None, '_aa4'),
                  on='aa4')

AAproperties.rename(columns={'aa4': 'aa5'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_merged, right=AAproperties, how='inner', indicator='location5', suffixes=(None, '_aa5'),
                  on='aa5')

AAproperties.rename(columns={'aa5': 'aa6'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_merged, right=AAproperties, how='inner', indicator='location6', suffixes=(None, '_aa6'),
                  on='aa6')

AAproperties.rename(columns={'aa6': 'aa7'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_merged, right=AAproperties, how='inner', indicator='location7', suffixes=(None, '_aa7'),
                  on='aa7')

AAproperties.rename(columns={'aa7': 'aa8'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_merged, right=AAproperties, how='inner', indicator='location8', suffixes=(None, '_aa8'),
                  on='aa8')

AAproperties.rename(columns={'aa8': 'aa9'}, inplace=True)
Ortho_merged = pd.merge(left=Ortho_merged, right=AAproperties, how='inner', indicator='location9', suffixes=(None, '_aa9'),
                  on='aa9')

X_test = Ortho_merged.drop(columns=['strain', 'locus', 'seq_type', 'aa_seq', 'pool_type', 'compound', 'median_s', 'Nham_aa', 'Resistance',
                                      'location1', 'location2', 'location3', 'location4', 'location5', 'location6', 'location7', 'location8', 'location9',
                                    'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9'])
y_test = Ortho_merged['Resistance']

########################################################################################################################

## Machine Learning

# Gridsearch & Random Forest
grid = {'n_estimators': [75, 100, 125, 150, 200, 250],
        'max_features': ['sqrt', 'log2', None],
        'max_depth': [5, 6, 7, None],
        'random_state': [18]
}

CV_rf = GridSearchCV(estimator=RandomForestClassifier(), param_grid=grid, n_jobs=-1, cv=10)
CV_rf.fit(X_train, y_train)

rf = CV_rf.best_estimator_
print(CV_rf.best_params_)

y_pred = rf.predict(X_test)

#rf = RandomForestClassifier(n_estimators= 250, max_features= 'sqrt', max_depth= 5, random_state= 18)
#rf.fit(X_train, y_train)
#y_pred = rf.predict(X_test)

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
plt.savefig(f'Conf_matrix.svg')
plt.show()

# shap
explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values[0], X_test, show=False, plot_size=(16,8), max_display=6)

plt.title(f'{drug} feature importance')
plt.gca().set_facecolor('white')
plt.grid(False)

# Modifying main plot parameters
plt.gca().tick_params(labelsize=20)
plt.gca().set_xlabel("SHAP value (impact on model output)", fontsize=20)
plt.gca().set_title(f'Caspofungin Feature Importance', fontsize=24)

# Get colorbar
cb_ax = plt.gcf().axes[1]

# Modifying color bar parameters
cb_ax.tick_params(labelsize=20)
cb_ax.set_ylabel("Feature value", fontsize=20)

# Adjust dot size
summary_plot = plt.gcf()
for scatter in summary_plot.axes[0].collections:
    scatter.set_sizes([100])  # Adjust the size as needed


# Define the custom y-axis tick labels you want to display for each feature
custom_ytick_labels = {
    'hydrophobicity_miyazawa': 'Position 1 Hydrophobicity Miyazawa',
    'hydrophobicity_guy': 'Position 1 Hydrophobicity Guy',
    'retention_coefficient_ph2.1_meek': 'Position 1 Retention Coefficient ph2.1 Meek',
    'percentage_accessible_residues_janin': 'Position 1 %Accessible Residues Janin',
    'hydrophobicity_antigenicity_welling_aa6': 'Position 6 Hydrophobicity Welling',
    'beta_sheet_levitt_aa9': 'Position 9 Beta Sheet Levitt',
}
# Create a function to get the custom label for a feature
def get_custom_label(feature):
    return custom_ytick_labels.get(feature, feature)

# Get the current y-axis tick labels
y_tick_labels = plt.gca().get_yticklabels()# Replace the y-axis tick labels with custom labels
custom_y_tick_labels = [get_custom_label(feature.get_text()) for feature in y_tick_labels]
plt.gca().set_yticklabels(custom_y_tick_labels)

plt.tight_layout()
plt.savefig(f'ShapValues.svg')
plt.show()

########################################################################################################################

# ROC curves
proba = rf.predict_proba(X_test)[::,1]
fpr, tpr, thresh = roc_curve(y_test, proba, pos_label='susceptible')
auc_score = roc_auc_score(y_test, proba)
print(f"AUC : {auc_score}")

random_probs = [0 for i in range(len(y_test))]
p_fpr, p_tpr, p_thresh = roc_curve(y_test, random_probs, pos_label='susceptible')

plt.figure(figsize=(15, 10))
plt.plot(fpr, tpr, color='blue', label= 'Random Forest Classifier')
plt.plot(p_fpr, p_tpr, linestyle='--', color='orange')
plt.title(f'Model ROC curve \n AUC : {auc_score:.3f}')
plt.xlabel('False Positive rate')
plt.ylabel('True positive rate')
plt.legend(loc= 'best')
# plt.tight_layout
plt.savefig('ROCcurve.svg')
plt.show()

"""