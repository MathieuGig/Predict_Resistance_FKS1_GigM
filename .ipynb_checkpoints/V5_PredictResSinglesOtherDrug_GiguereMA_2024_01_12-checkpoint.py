## Title: PredictResSinglesOtherDrug
## Author: Mathieu Giguere
## Brief: Uses machine learning to predict the resistance of FKS1-HS1 single mutants amino acid sequences to an
##        antifungal drug by training a model on the same single mutants amino acid sequences exposed to another drug.
##        The model uses Expasy Protscale's amino acid properties as features.
## Preconditions: Needs Romain Durand's 'DMS-main' repository and 'aminoAcidProperties.txt' and 'BY11nnkNT_classes.csv'.

# importing modules and packages.
import pandas as pd
import numpy as np
np.bool = np.bool_
np.int = np.int_
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, matthews_corrcoef, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
import shap  # v0.39.0
shap.initjs()

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

########################################################################################################################

# VIP. Very important parameter. The antifungal drugs.
#drug1 = 'caspofungin'
#drug1 = 'micafungin'
drug1 = 'anidulafungin'

#drug2 = 'caspofungin'
drug2 = 'micafungin'
#drug2 = 'anidulafungin'

########################################################################################################################

# importing data from Romain Durand's 'BY11nnkNT_classes.csv'
df = pd.read_csv('BY11nnkNT_classes.csv')

VerifyConditions_df = df.loc[df['compound'].isin(['anidulafungin', 'caspofungin', 'micafungin']) &
                             (df['seq_type'] == 'single') &
                             (df['pool_type'] == 'single')].reset_index()

# Drop some columns to simplify the dataframe.
simplified_df = VerifyConditions_df.drop(columns=['index', 'Unnamed: 0', 'seq_type', 'pool_type', 'Nham_aa', 'cscore'])

# intermediary is considered resistant.
simplified_df['class'] = np.where((simplified_df['class'] == 'intermediary'), 'resistant', simplified_df['class'])

########################################################################################################################

# Amino acid properties
AAproperties = pd.read_table('aminoAcidProperties.txt')
AAproperties.rename(columns={'Aminoacid.1.letter': 'aa0'}, inplace=True)

########################################################################################################################

# Single mutations dataframe on drug 1.
Drug1_df = simplified_df.loc[simplified_df['compound'] == drug1]

# Get a column for each amino acid of the aa_seq.
Drug1_df[['aa0', 'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8']] = Drug1_df['aa_seq'].apply(
    lambda x: pd.Series(list(x)))

# Get rid of '*' mutants.
Drug1_df = Drug1_df.loc[(Drug1_df['aa0'] != '*') & (Drug1_df['aa1'] != '*') & (Drug1_df['aa2'] != '*') &
                        (Drug1_df['aa3'] != '*') & (Drug1_df['aa4'] != '*') & (Drug1_df['aa5'] != '*') &
                        (Drug1_df['aa6'] != '*') & (Drug1_df['aa7'] != '*') & (Drug1_df['aa8'] != '*')]

# Perform a series of merge operations to associate the amino acid properties to the corresponding sequence.
Drug1_merged = pd.merge(left=Drug1_df, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa0'),
                        on='aa0')

AAproperties.rename(columns={'aa0': 'aa1'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa1'),
                        on='aa1')

AAproperties.rename(columns={'aa1': 'aa2'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa2'),
                        on='aa2')

AAproperties.rename(columns={'aa2': 'aa3'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa3'),
                        on='aa3')

AAproperties.rename(columns={'aa3': 'aa4'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa4'),
                        on='aa4')

AAproperties.rename(columns={'aa4': 'aa5'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa5'),
                        on='aa5')

AAproperties.rename(columns={'aa5': 'aa6'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa6'),
                        on='aa6')

AAproperties.rename(columns={'aa6': 'aa7'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa7'),
                        on='aa7')

AAproperties.rename(columns={'aa7': 'aa8'}, inplace=True)
Drug1_merged = pd.merge(left=Drug1_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa8'),
                        on='aa8')

# Get training data for machine learning.
X_train = Drug1_merged.drop(columns=['compound', 'aa_seq', 'nt_seq', 's_corr', 'class',
                                     'aa0', 'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8'])

y_train = Drug1_merged['class']

########################################################################################################################

# Single mutations dataframe on drug 2
Drug2_df = simplified_df.loc[simplified_df['compound'] == drug2]

# Get a column for each amino acid of the aa_seq.
Drug2_df[['aa0', 'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8']] = Drug2_df['aa_seq'].apply(
    lambda x: pd.Series(list(x)))

# Get rid of '*' mutants.
Drug2_df = Drug2_df.loc[(Drug2_df['aa0'] != '*') & (Drug2_df['aa1'] != '*') & (Drug2_df['aa2'] != '*') &
                        (Drug2_df['aa3'] != '*') & (Drug2_df['aa4'] != '*') & (Drug2_df['aa5'] != '*') &
                        (Drug2_df['aa6'] != '*') & (Drug2_df['aa7'] != '*') & (Drug2_df['aa8'] != '*')]

# Perform a series of merge operations to associate the amino acid properties to the corresponding sequence.
AAproperties.rename(columns={'aa8': 'aa0'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_df, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa0'),
                        on='aa0')

AAproperties.rename(columns={'aa0': 'aa1'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa1'),
                        on='aa1')

AAproperties.rename(columns={'aa1': 'aa2'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa2'),
                        on='aa2')

AAproperties.rename(columns={'aa2': 'aa3'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa3'),
                        on='aa3')

AAproperties.rename(columns={'aa3': 'aa4'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa4'),
                        on='aa4')

AAproperties.rename(columns={'aa4': 'aa5'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa5'),
                        on='aa5')

AAproperties.rename(columns={'aa5': 'aa6'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa6'),
                        on='aa6')

AAproperties.rename(columns={'aa6': 'aa7'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa7'),
                        on='aa7')

AAproperties.rename(columns={'aa7': 'aa8'}, inplace=True)
Drug2_merged = pd.merge(left=Drug2_merged, right=AAproperties, how='inner', indicator=False, suffixes=(None, '_aa8'),
                        on='aa8')

# Get test data for machine learning.
X_test = Drug2_merged.drop(columns=['compound', 'aa_seq', 'nt_seq', 's_corr', 'class',
                                    'aa0', 'aa1', 'aa2', 'aa3', 'aa4', 'aa5', 'aa6', 'aa7', 'aa8'])

y_test = Drug2_merged['class']

########################################################################################################################
# Machine Learning
# Use Gridsearch & Random Forest.
grid = {'n_estimators': [150, 200, 250, 300, 400, 500],
        'max_features': ['sqrt', 'log2', None],
        'max_depth': [5, 6, 7, 10],
        'random_state': [18]}
CV_rf = GridSearchCV(estimator=RandomForestClassifier(), param_grid=grid, n_jobs=-1, cv=10)
CV_rf.fit(X_train, y_train)
rf = CV_rf.best_estimator_
print(CV_rf.best_params_)

rf.fit(X_train, y_train)

y_pred = rf.predict(X_test)

# Calculate accuracy of predictions.
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)

# Calculate Matthews correlation coefficient.
mat = matthews_corrcoef(y_test, y_pred)
print("Mat:", mat)

# Create the confusion matrix.
cmatrix = confusion_matrix(y_test, y_pred, labels=['resistant', 'WT-like'])

# Display the confusion matrix using seaborn heatmap.
plt.figure(figsize=(10, 10))
sns.set(font_scale=3)
sns.heatmap(cmatrix, annot=True, fmt='d', cmap='Blues', xticklabels=['resistant', 'WT-like'],
            yticklabels=['resistant', 'WT-like'], annot_kws={'size': 60})
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title(f'Train: {drug1} \n Pred: {drug2} \n Confusion Matrix \n Accuracy : {accuracy:.3f}')
plt.tight_layout()
plt.savefig(f'V4_ConfMatrix_{drug1}_toPredict_{drug2}.svg')
plt.show()

########################################################################################################################

# Perform shap values analysis and display plot.
explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values[1], X_test, show=False, plot_size=(16, 8),
                  max_display=6)  # shap values[1] targets resistance

# Modify background of plot.
plt.gca().set_facecolor('white')
plt.grid(False)

# Modifying main plot parameters.
plt.gca().tick_params(labelsize=20)
plt.gca().set_xlabel("SHAP value (impact on model output)", fontsize=20)
plt.gca().set_title(f'Train: {drug1} \n Pred: {drug2} \n Feature Importance', fontsize=24)

# Get colorbar & modify color bar parameters.
cb_ax = plt.gcf().axes[1]
cb_ax.tick_params(labelsize=20)
cb_ax.set_ylabel("Feature value", fontsize=20)

# Adjust dot size.
summary_plot = plt.gcf()
for scatter in summary_plot.axes[0].collections:
    scatter.set_sizes([100])  # Adjust the size as needed

plt.tight_layout()
plt.savefig(f'V4_ShapValues_TrainOn_{drug1}_toPredict_{drug2}.svg')
plt.show()

########################################################################################################################

# Create AUC-ROC curves.
proba = rf.predict_proba(X_test)[::, 1]
fpr, tpr, thresh = roc_curve(y_test, proba, pos_label='resistant')
auc_score = roc_auc_score(y_test, proba)

random_probs = [0 for i in range(len(y_test))]
p_fpr, p_tpr, p_thresh = roc_curve(y_test, random_probs, pos_label='resistant')

# Make plot.
plt.figure(figsize=(15, 15))
plt.plot(fpr, tpr, color='blue', label='Random Forest Classifier')
plt.plot(p_fpr, p_tpr, linestyle='--', color='orange')
plt.title(f'Train: {drug1}. Pred: {drug2} \n Model ROC curve \n AUC : {auc_score:.3f}')
plt.xlabel('False Positive rate')
plt.ylabel('True positive rate')
plt.legend(loc='best')
plt.savefig(f'V4_ROC_TrainOn_{drug1}_toPredict_{drug2}.svg')
plt.show()

########################################################################################################################
# Heatmaps

SummaryData = Drug2_merged[["aa_seq", "s_corr", "class"]]

# Find misclassified values.
misclassified = np.where(y_pred != y_test)
List_misclassified = misclassified[0].tolist()
SummaryData['misclassified'] = np.where(SummaryData.index.isin(List_misclassified), 1, None)

# recreate 'aa_pos' and 'alt_aa' columns to make heatmap.
# where 'aa_pos' is the position of the mutant amino acid in the hotspot sequence.
# and where 'alt_aa' is the one letter code for the mutant amino acid.

def difference(string1, string2):
    # this function compares two strings and returns the characters unique to string1.
    # Unpack strings into lists.
    string1 = [*string1]
    string2 = [*string2]

    setA = set(string1)  # Store all string list items in set A.
    setB = set(string2)

    str_diff = setA.difference(setB)  # returns items that are unique to the first set.
    isEmpty = (len(str_diff) == 0)

    if isEmpty:
        return 'n'
    else:
        return str_diff

# Initilise lists and transform WT sequence into list of amino acids.
alt_aa_list = []
pos_aa_list = []
wtaa = 'FLVLSLRDP'  # FKS1-HS1 Hotspot in aa
split_wt = [*wtaa]

# Compare each sequence to the WT and find the position and amino acid of each mutant.
# Append None if sequence is not a mutant.
for val in SummaryData['aa_seq']:
    split_seq = [*val]
    is_mutant = False
    for num in range(9):
        dif = list(difference(split_seq[num], split_wt[num]))[0]
        if dif != 'n':
            alt_aa_list.append(dif)
            pos_aa_list.append(num)
            is_mutant = True
            break
    if not is_mutant:
        alt_aa_list.append(None)
        pos_aa_list.append(None)

# Add 'alt_aa' and 'aa_pos' columns to dataframe for heatmaps.
SummaryData['alt_aa'] = alt_aa_list
SummaryData['aa_pos'] = pos_aa_list

# Construct heatmap dataframe.
heatmapDF_single = SummaryData.groupby(['alt_aa', 'aa_pos'])[['s_corr']].median().reset_index()
heatmapDF_wide = heatmapDF_single.pivot(index='alt_aa', columns='aa_pos', values='s_corr')

heatmapDF_mask = SummaryData.groupby(['alt_aa', 'aa_pos'])[['misclassified']].max().reset_index()
heatmapDF_wide_mask = heatmapDF_mask.pivot(index='alt_aa', columns='aa_pos', values='misclassified')

# Sort dataframe to match Romain Durand's.
aa_sort_order = 'PGCQNTSEDKHRWYFMLIVA'
aa_sort_dic = dict(zip(list(aa_sort_order), list(range(0, len(aa_sort_order)))))
heatmapDF_wide.sort_index(key=lambda x: x.map(aa_sort_dic), inplace=True)
heatmapDF_wide_mask.sort_index(key=lambda x: x.map(aa_sort_dic), inplace=True)

# Purge graph space.
sns.set(rc={'figure.figsize': (4, 6), })
sns.set_theme(style='white')
f, ax = plt.subplots()

# Romain Durand's custom color palette.
ccmap = sns.color_palette("blend:#009B9E,#42B7B9,#A7D3D4,#F1F1F1,#E4C1D9,#D691C1,#C75DAB",  # CARTOColors Tropic
                          as_cmap=True)
ccmap.set_bad('.5')  # Color for missing values or masked.

# Mask the correctly classified values.
mask = pd.isnull(heatmapDF_wide_mask)

# Draw heatmap.
ax = sns.heatmap(heatmapDF_wide, mask=mask,
                 cmap=ccmap, center=0,
                 vmin=-2,
                 vmax=2)
ax.set_xlabel('Position in hotspot')
ax.set_ylabel(None)
plt.yticks(rotation=0)
ax.set_title(f'Train on {drug1} to predict {drug2}.\n Heatmap of misclassified mutants')

# Here, missing values from the data are put in black.
for i, j in zip(*np.where(heatmapDF_wide.isnull())):
    plt.gca().add_patch(plt.Rectangle((j, i), 1, 1, fill=True, color='black'))

wtaa = 'FLVLSLRDP'  # FKS1-HS1 Hotspot in aa
wtcoord_aa = [(i + 0.5, list(aa_sort_order).index(v) + 0.5) for i, v in enumerate(wtaa)]

for o in wtcoord_aa:
    ax.plot(o[0], o[1], 'or')  # displays WT codons by marking them with red (r) circles (o)

plt.savefig(f'V4_Heatmap_misclassified_TrainOn_{drug1}_ToPredict_{drug2}.svg', format='svg', dpi=300)
plt.show()
