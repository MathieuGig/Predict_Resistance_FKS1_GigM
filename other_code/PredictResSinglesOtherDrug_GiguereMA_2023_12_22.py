## Title: PredictResSinglesOtherDrug
## Author: Mathieu Giguere
## Brief: Uses machine learning to predict the resistance of FKS1-HS1 single mutants amino acid sequences to an
##        antifungal drug by training a model on the same single mutants amino acid sequences exposed to another drug.
##        The model uses Expasy Protscale's amino acid properties as features.
## Preconditions: Needs Romain Durand's 'DMS-main' repository and 'aminoAcidProperties.txt'.

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
#drug1 = 'caspofungin'
#drug1 = 'micafungin'
drug1 = 'anidulafungin'

drug2 = 'caspofungin'
#drug2 = 'micafungin'
#drug2 = 'anidulafungin'

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
AAproperties = pd.read_table('aminoAcidProperties.txt')
AAproperties.rename(columns={'Aminoacid.1.letter': 'aa1'}, inplace=True)

########################################################################################################################

# Single mutations dataframe on drug 1
Drug1_master = shortMaster.loc[(shortMaster['seq_type'] == 'single') & (shortMaster['pool_type'] == 'single') & (shortMaster['strain'] == 'BY4741') & (shortMaster['locus'] == 'FKS1-HS1') & (shortMaster['compound'] == drug1)]
Drug1_master = Drug1_master.loc[(Drug1_master['Nham_aa'] == 1) | ((Drug1_master['Nham_aa'] == 0) & (Drug1_master['median_s'] < 1))]
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
Drug2_master = Drug2_master.loc[(Drug2_master['Nham_aa'] == 1) | ((Drug2_master['Nham_aa'] == 0) & (Drug2_master['median_s'] < 1))]
Drug2_master['Resistance'] = np.where(Drug2_master['median_s'] >= 1, 'resistant', 'susceptible')

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

X_test = Drug2_merged.drop(columns=['strain', 'locus', 'seq_type', 'aa_seq', 'pool_type', 'compound', 'aa_pos', 'alt_aa', 'median_s', 'Nham_aa', 'Resistance',
                                      'location1', 'location2', 'location3', 'location4', 'location5', 'location6', 'location7', 'location8', 'location9',
                                    'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9'])
y_test = Drug2_merged['Resistance']

########################################################################################################################
## Machine Learning

# Gridsearch & Random Forest
#grid = {'n_estimators': [200, 300, 400, 500, 600],
#        'max_features': ['sqrt', 'log2', None],
#       'max_depth': [5, 6, 7, None],
#        'random_state': [18]
#}
#CV_rf = GridSearchCV(estimator=RandomForestClassifier(), param_grid=grid, n_jobs=-1, cv=10)
#CV_rf.fit(X_train, y_train)
#rf = CV_rf.best_estimator_
#print(CV_rf.best_params_)

# rf best parameters by drugs according to GridSearchCV.
# Train: caspofungin. Pred: anidulafungin
# Train: caspofungin. Pred: micafungin
# Train: anidulafungin. Pred: caspofungin
# Train: anidulafungin. Pred: micafungin
rf = RandomForestClassifier(n_estimators= 500, max_features= 'sqrt', max_depth= 5, random_state= 18)

# Train: micafungin. Pred: caspofungin
# Train: micafungin. Pred: anidulafungin
#rf = RandomForestClassifier(n_estimators= 200, max_features= 'sqrt', max_depth= None, random_state= 18)

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
plt.title(f'Train: {drug1} \n Pred: {drug2} \n Confusion Matrix \n Accuracy : {accuracy:.3f}')
plt.tight_layout()
plt.savefig(f'ConfMatrix_TrainOn_{drug1}_toPredict_{drug2}.svg')
plt.show()

########################################################################################################################

# shap
explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values[0], X_test, show=False, plot_size=(16,8), max_display=6)

# Modify background
plt.gca().set_facecolor('white')
plt.grid(False)

# Modifying main plot parameters
plt.gca().tick_params(labelsize=20)
plt.gca().set_xlabel("SHAP value (impact on model output)", fontsize=20)
plt.gca().set_title(f'Train: {drug1} \n Pred: {drug2} \n Feature Importance', fontsize=24)

# Get colorbar
cb_ax = plt.gcf().axes[1]

# Modifying color bar parameters
cb_ax.tick_params(labelsize=20)
cb_ax.set_ylabel("Feature value", fontsize=20)

# Adjust dot size
#summary_plot = plt.gcf()
#for scatter in summary_plot.axes[0].collections:
#    scatter.set_sizes([100])  # Adjust the size as needed

plt.tight_layout()
plt.savefig(f'ShapValues_TrainOn_{drug1}_toPredict_{drug2}.svg')
plt.show()


########################################################################################################################

# ROC curves
proba = rf.predict_proba(X_test)[::,1]
fpr, tpr, thresh = roc_curve(y_test, proba, pos_label='susceptible')
auc_score = roc_auc_score(y_test, proba)
#print(f"AUC : {auc_score}")

random_probs = [0 for i in range(len(y_test))]
p_fpr, p_tpr, p_thresh = roc_curve(y_test, random_probs, pos_label='susceptible')

plt.figure(figsize=(15, 15))
plt.plot(fpr, tpr, color='blue', label= 'Random Forest Classifier')
plt.plot(p_fpr, p_tpr, linestyle='--', color='orange')
plt.title(f'Train: {drug1}. Pred: {drug2} \n Model ROC curve \n AUC : {auc_score:.3f}')
plt.xlabel('False Positive rate')
plt.ylabel('True positive rate')
plt.legend(loc='best')
plt.savefig(f'ROC_TrainOn_{drug1}_toPredict_{drug2}.svg')
plt.show()

########################################################################################################################

# Heatmaps
# dataframe avec alt_aa en rows, aa_pos en colonnes, puis values == misclassified ?
SummaryData = Drug2_merged[["aa_seq", "median_s", "alt_aa", "aa_pos" ,"Resistance"]]

# Find misclassified values
misclassified = np.where(y_pred != y_test)
List_misclassified = misclassified[0].tolist()
SummaryData['misclassified'] = np.where(SummaryData.index.isin(List_misclassified), 1, None)

# Construct heatmap dataframe
heatmapDF_single = SummaryData.groupby(['alt_aa','aa_pos'])[['median_s']].median().reset_index()
heatmapDF_wide = heatmapDF_single.pivot(index='alt_aa', columns='aa_pos', values='median_s')

heatmapDF_mask = SummaryData.groupby(['alt_aa','aa_pos'])[['misclassified']].max().reset_index()
heatmapDF_wide_mask = heatmapDF_mask.pivot(index='alt_aa', columns='aa_pos', values='misclassified')

# Purge graph space
sns.set(rc = {'figure.figsize':(4,6),
             })
sns.set_theme(style='white')
f, ax = plt.subplots()

# Custom color palette
ccmap = sns.color_palette("blend:#009B9E,#42B7B9,#A7D3D4,#F1F1F1,#E4C1D9,#D691C1,#C75DAB", # CARTOColors Tropic
                          as_cmap=True)
ccmap.set_bad('.5') # Color for missing values or masked

# Mask the correctly classified values.
mask = pd.isnull(heatmapDF_wide_mask)

# Draw heatmap
ax = sns.heatmap(heatmapDF_wide, mask=mask,
                 cmap=ccmap, center = 0,
                 vmin=-2,
                 vmax=2
                )
ax.set_xlabel('Position in hotspot')
ax.set_ylabel(None)
plt.yticks(rotation=0)
ax.set_title(f'Train on {drug1} to predict {drug2}.\n Heatmap of misclassified mutants')

aa_sort_order = 'ACEDFGHIKLMNPQRSTVWY' # alphabetical
wtaa = 'FLVLSLRDP' #FKS1-HS1 Hotspot in aa
wtcoord_aa = [(i+0.5, list(aa_sort_order).index(v)+0.5) for i,v in enumerate(wtaa)]

for o in wtcoord_aa:
    ax.plot(o[0],o[1],'or') # displays WT codons by marking them with red (r) circles (o)

plt.savefig(f'Heatmap_misclassified_TrainOn_{drug1}_ToPredict_{drug2}.svg', format='svg', dpi=300)
plt.show()