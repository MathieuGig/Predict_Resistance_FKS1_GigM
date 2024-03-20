# Predict_Resistance_FKS1
Mathieu Giguère's Summer 2023 internship at the LandryLab and his work during the Fall 2023 semester.\
Most files require to be in the same directory as Romain Durand's 'DMS-main' repository. The code in my repository is made to work with it. Be careful with absolute and relative file paths.

Folders :
- 'other_code' : All previous saves of my code. Mostly uncommented and unattended to.

Files :
- 'aminoAcidProperties.txt' : The amino acid properties used as features for machine learning. They come from Expasy's Protscale tool. \
Link: https://web.expasy.org/protscale/
- 'PredictResistance_GiguereMA_2023_11_29.py' : Python code to predict the resistance of FKS1 orthologs to echinocandins.
- 'V5_PredictResSinglesOtherDrug_GiguereMA_2024_01_12.py' : Python code to predict the resistance of FKS1 single mutants to different echinocandins.
- 'fks1_dataframe.csv' : A result of my previous ortholog taxonomic analysis. I use it in 'PredictResistance_GiguereMA_2023_11_29.py'.
- 'BY11nnkNT_classes.csv' : Romain Durand's data on fks1 single mutants. I use it in 'V5_PredictResSinglesOtherDrug_GiguereMA_2024_01_12.py'.

- 'Predict_orthologs.ipynb' : jupyter notebook to predict the resistance of FKS1 orthologs to echinocandins.
- 'Make_cmatrix.ipnyb' : jupyter notebook to make a confusion matrix based on the output of 'Predict_orthologs.ipynb'.
- 'Make_ROCcurve.ipnyb' : jupyter notebook to make a ROC curve based on the output of 'Predict_orthologs.ipynb'.
