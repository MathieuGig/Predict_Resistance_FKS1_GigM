import pandas as pd
import numpy as np
import re

from GenerateWebLogo import GenerateWebLogo

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

def makeLogoPairPlots(strain, locus):
    df = pd.read_csv(f'Data_pairPlot_{strain}_{locus}_single_ortho_.csv')
    df['Resistant_Anidulafungin'] = np.where(df['Anidulafungin'] >= 1.5, True, False)
    df['Resistant_Caspofungin'] = np.where(df['Caspofungin'] >= 1.5, True, False)
    df['Resistant_Micafungin'] = np.where(df['Micafungin'] >= 1.5, True, False)

    R_ALL = df.loc[(df['Resistant_Anidulafungin'] == True) & (df['Resistant_Caspofungin'] == True) & (df['Resistant_Micafungin'] == True)]
    R_ALL_sequences = np.array2string(R_ALL['Sequence'].values)

    R_Ani_Caspo = df.loc[(df['Resistant_Anidulafungin'] == True) & (df['Resistant_Caspofungin'] == True) & (df['Resistant_Micafungin'] == False)]
    R_Ani_Caspo_sequences = np.array2string(R_Ani_Caspo['Sequence'].values)

    R_Ani_Mica = df.loc[(df['Resistant_Anidulafungin'] == True) & (df['Resistant_Caspofungin'] == False) & (df['Resistant_Micafungin'] == True)]
    R_Ani_Mica_sequences = np.array2string(R_Ani_Mica['Sequence'].values)

    R_Caspo_Mica = df.loc[(df['Resistant_Anidulafungin'] == False) & (df['Resistant_Caspofungin'] == True) & (df['Resistant_Micafungin'] == True)]
    R_Caspo_Mica_sequences = np.array2string(R_Caspo_Mica['Sequence'].values)

    R_OnlyAni = df.loc[(df['Resistant_Anidulafungin'] == True) & (df['Resistant_Caspofungin'] == False) & (df['Resistant_Micafungin'] == False)]
    R_OnlyAni_sequences = np.array2string(R_OnlyAni['Sequence'].values)

    R_OnlyCaspo = df.loc[(df['Resistant_Anidulafungin'] == False) & (df['Resistant_Caspofungin'] == True) & (df['Resistant_Micafungin'] == False)]
    R_OnlyCaspo_sequences = np.array2string(R_OnlyCaspo['Sequence'].values)

    R_OnlyMica = df.loc[(df['Resistant_Anidulafungin'] == False) & (df['Resistant_Caspofungin'] == False) & (df['Resistant_Micafungin'] == True)]
    R_OnlyMica_sequences = np.array2string(R_OnlyMica['Sequence'].values)

    if len(R_ALL_sequences) > 13:
        GenerateWebLogo(R_ALL_sequences, f'{strain}_{locus}_all_resistant')

    if len(R_Ani_Caspo_sequences) > 13:
        GenerateWebLogo(R_Ani_Caspo_sequences, f'{strain}_{locus}_Ani_Caspo_resistant')

    if len(R_Ani_Mica_sequences) > 13:
        GenerateWebLogo(R_Ani_Mica_sequences, f'{strain}_{locus}_Ani_Mica_resistant')

    if len(R_Caspo_Mica_sequences) > 13:
        GenerateWebLogo(R_Caspo_Mica_sequences, f'{strain}_{locus}_Caspo_Mica_resistant')

    if len(R_OnlyAni_sequences) > 13:
        GenerateWebLogo(R_OnlyAni_sequences, f'{strain}_{locus}_Ani_resistant')

    if len(R_OnlyCaspo_sequences) > 13:
        GenerateWebLogo(R_OnlyCaspo_sequences, f'{strain}_{locus}_Caspo_resistant')

    if len(R_OnlyMica_sequences) > 13:
        GenerateWebLogo(R_OnlyMica_sequences, f'{strain}_{locus}_Mica_resistant')

makeLogoPairPlots('BY4741', 'FKS1-HS1')
makeLogoPairPlots('BY4741', 'FKS1-HS2')
makeLogoPairPlots('R1158', 'FKS1-HS1')
makeLogoPairPlots('R1158', 'FKS1-HS2')
makeLogoPairPlots('R1158', 'FKS2-HS1')
makeLogoPairPlots('R1158', 'FKS2-HS2')