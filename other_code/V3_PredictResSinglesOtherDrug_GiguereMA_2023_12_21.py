## Title: PredictResSinglesOtherDrug
## Author: Mathieu Giguere
## Brief: Uses machine learning to predict the resistance of FKS1-HS1 single mutants amino acid sequences to an
##        antifungal drug by training a model on the same single mutants amino acid sequences exposed to another drug.
##        The model uses Expasy Protscale's amino acid properties as features.
## Preconditions: Needs Romain Durand's 'DMS-main' repository and 'aminoAcidProperties.txt'. and
##                the 'GenerateWebLogo' function from the 'GenerateWebLogo_GiguereMA_2023_11_29'
##                file.

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

# importing data
df = pd.read_csv('BY11nnkNT_classes.csv')

VerifyConditions_df = df.loc[df['compound'].isin(['anidulafungin', 'caspofungin', 'micafungin']) &
                             (df['seq_type'] == 'single') &
                             (df['pool_type'] == 'single')].reset_index()

simplified_df = VerifyConditions_df.drop(columns=['index', 'Unnamed: 0', 'seq_type', 'pool_type', 'Nham_aa', 'cscore'])

########################################################################################################################

# Amino acid properties
AAproperties = pd.read_table('aminoAcidProperties.txt')
AAproperties.rename(columns={'Aminoacid.1.letter': 'aa1'}, inplace=True)

########################################################################################################################

# Single mutations dataframe on drug 1
Drug1_df = simplified_df.loc[simplified_df['compound'] == drug1]

Drug1_df[['aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9']] = Drug1_df['aa_seq'].apply(lambda x: pd.Series(list(x)))
Drug1_df = Drug1_df.loc[(Drug1_df['aa1'] != '*') & (Drug1_df['aa2'] != '*') & (Drug1_df['aa3'] != '*') &
                        (Drug1_df['aa4'] != '*') & (Drug1_df['aa5'] != '*') & (Drug1_df['aa6'] != '*') &
                        (Drug1_df['aa7'] != '*') & (Drug1_df['aa8'] != '*') & (Drug1_df['aa9'] != '*')]

Drug1_merged = pd.merge(left=Drug1_df, right=AAproperties, how='inner', indicator='location1', suffixes=(None, '_aa1'),
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

X_train = Drug1_merged.drop(columns=['compound', 'aa_seq', 'nt_seq', 's_corr', 'class',
                                     'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9',
                                     'location1', 'location2', 'location3', 'location4', 'location5',
                                     'location6', 'location7', 'location8', 'location9',])

y_train = Drug1_merged['class']

########################################################################################################################

# Single mutations dataframe on drug 2
Drug2_df = simplified_df.loc[simplified_df['compound'] == drug2]

Drug2_df[['aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9']] = Drug2_df['aa_seq'].apply(lambda x: pd.Series(list(x)))
Drug2_df = Drug2_df.loc[(Drug2_df['aa1'] != '*') & (Drug2_df['aa2'] != '*') & (Drug2_df['aa3'] != '*') &
                        (Drug2_df['aa4'] != '*') & (Drug2_df['aa5'] != '*') & (Drug2_df['aa6'] != '*') &
                        (Drug2_df['aa7'] != '*') & (Drug2_df['aa8'] != '*') & (Drug2_df['aa9'] != '*')]

AAproperties.rename(columns={'aa9': 'aa1'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_df, right=AAproperties, how='inner', indicator='location1', suffixes=(None, '_aa1'),
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

X_test = Drug2_merged.drop(columns=['compound', 'aa_seq', 'nt_seq', 's_corr', 'class',
                                     'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8', 'aa9',
                                    'location1', 'location2', 'location3', 'location4', 'location5',
                                     'location6', 'location7', 'location8', 'location9'])

y_test = Drug2_merged['class']

########################################################################################################################
## Machine Learning
#print('go ML')
# Gridsearch & Random Forest
#grid = {'n_estimators': [100, 125, 150, 200, 250, 300, 400, 500, 600],
#        'max_features': ['sqrt', 'log2', None],
#        'max_depth': [5, 6, 7, 10],
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
rf = RandomForestClassifier(n_estimators= 400, max_features= 'sqrt', max_depth= 10, random_state= 18)

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


# Create the confusion matrix
cmatrix = confusion_matrix(y_test, y_pred)

# Display the confusion matrix using seaborn heatmap
plt.figure(figsize=(10, 10))
sns.set(font_scale=3)
sns.heatmap(cmatrix, annot=True, fmt='d', cmap='Blues', xticklabels=['resistant', 'intermediary', "WT-like"], yticklabels=['resistant', 'intermediary', "WT-like"], annot_kws={"size": 60})
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title(f'Train: {drug1} \n Pred: {drug2} \n Confusion Matrix \n Accuracy : {accuracy:.3f}')
plt.tight_layout()
plt.savefig(f'ConfMatrix_{drug1}_toPredict_{drug2}.svg')
plt.show()

########################################################################################################################
"""
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

## WebLogo generation of the confusion Matrix sequences.

# Make dataframe with these columns: sequence, true label, classification label, misclassifed.

SummaryData = Drug2_merged[["aa_seq", "Resistance"]]

# Find misclassified
misclassified = np.where(y_pred != y_test)
List_misclassified = misclassified[0].tolist()

SummaryData['misclassified'] = np.where(SummaryData.index.isin(List_misclassified), True, False)
SummaryData['Prediction'] = np.where((SummaryData['misclassified'] == True) &
                                     (SummaryData['Resistance'] == 'susceptible'), 'resistant',
                                     np.where((SummaryData['misclassified'] == True) &
                                              (SummaryData['Resistance'] == 'resistant'),
                                              'susceptible', SummaryData['Resistance']))

# Use code from GenerateWebLogos
# Define conditions
# R = Resistant. S = Susceptible. RR = sequence resistant et prediction resistante
# Put sequences into string
# Call GenerateWebLogo
# It prints a command you can use in a terminal to generate WebLogos using 'WebLogo 3'.
df_RR = SummaryData.loc[(SummaryData['Resistance'] == 'resistant') & (SummaryData['Prediction'] == 'resistant')]
Seq_RR = np.array2string(df_RR['aa_seq'].values)
GenerateWebLogo(Seq_RR, f'RR_{drug1}_{drug2}')

df_SS = SummaryData.loc[(SummaryData['Resistance'] == 'susceptible') & (SummaryData['Prediction'] == 'susceptible')]
Seq_SS = np.array2string(df_SS['aa_seq'].values)
GenerateWebLogo(Seq_SS, f'SS_{drug1}_{drug2}')

df_RS = SummaryData.loc[(SummaryData['Resistance'] == 'resistant') & (SummaryData['Prediction'] == 'susceptible')]
Seq_RS = np.array2string(df_RS['aa_seq'].values)
GenerateWebLogo(Seq_RS, f'RS_{drug1}_{drug2}')

df_SR = SummaryData.loc[(SummaryData['Resistance'] == 'susceptible') & (SummaryData['Prediction'] == 'resistant')]
Seq_SR = np.array2string(df_SR['aa_seq'].values)
GenerateWebLogo(Seq_SR, f'SR_{drug1}_{drug2}')

print('\n Copy paste the commands into a terminal to use WebLogo 3 and generate the sequence logos.')


########################################################################################################################

# Heatmaps
# dataframe avec alt_aa en rows, aa_pos en colonnes, puis values == misclassified ?
# normaliser sur le wildtype
# Make dataframe puis plug dans seaborn heatmap
SummaryData = Drug2_merged[["aa_seq", "alt_aa", "aa_pos" ,"Resistance"]]

# Find misclassified
misclassified = np.where(y_pred != y_test)
List_misclassified = misclassified[0].tolist()

SummaryData['misclassified'] = np.where(SummaryData.index.isin(List_misclassified), 1, 0)
SummaryData['Prediction'] = np.where((SummaryData['misclassified'] == 1) &
                                     (SummaryData['Resistance'] == 'susceptible'), 'resistant',
                                     np.where((SummaryData['misclassified'] == 1) &
                                              (SummaryData['Resistance'] == 'resistant'),
                                              'susceptible', SummaryData['Resistance']))

heatmapDF_single = SummaryData.groupby(['alt_aa','aa_pos'])[['misclassified']].median().reset_index()
heatmapDF_wide = heatmapDF_single.pivot(index='alt_aa', columns='aa_pos', values='misclassified')

# Purge graph space
sns.set(rc = {'figure.figsize':(4,6),
             })
sns.set_theme(style='white')
f, ax = plt.subplots()

# Custom color palette
ccmap = sns.color_palette("blend:#009B9E,#42B7B9,#A7D3D4,#F1F1F1,#E4C1D9,#D691C1,#C75DAB", # CARTOColors Tropic
                          as_cmap=True)
ccmap.set_bad('.5') # Color for missing values

# Mask for missing values
mask = pd.isnull(heatmapDF_wide)

# Draw heatmap
ax = sns.heatmap(heatmapDF_wide, mask=mask,
                 cmap=ccmap, center = 0.5,
                 vmin=0,
                 vmax=1
                )
ax.set_xlabel('Position in hotspot')
ax.set_ylabel(None)
plt.yticks(rotation=0)
ax.set_title('Nice')

plt.savefig('aa_heatmap_misclassified.svg', format='svg', dpi=300)
"""