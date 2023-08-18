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

    V_ALL = df.loc[(df['Resistant_Anidulafungin'] == False) & (df['Resistant_Caspofungin'] == False) & (df['Resistant_Micafungin'] == False)]
    V_ALL_sequences = np.array2string(V_ALL['Sequence'].values)

    R_Ani_Caspo = df.loc[(df['Resistant_Anidulafungin'] == True) & (df['Resistant_Caspofungin'] == True)]
    R_Ani_Caspo_sequences = np.array2string(R_Ani_Caspo['Sequence'].values)

    R_Ani_Mica = df.loc[(df['Resistant_Anidulafungin'] == True) & (df['Resistant_Micafungin'] == True)]
    R_Ani_Mica_sequences = np.array2string(R_Ani_Mica['Sequence'].values)

    R_Caspo_Mica = df.loc[(df['Resistant_Caspofungin'] == True) & (df['Resistant_Micafungin'] == True)]
    R_Caspo_Mica_sequences = np.array2string(R_Caspo_Mica['Sequence'].values)


    V_Ani_Caspo = df.loc[(df['Resistant_Anidulafungin'] == False) & (df['Resistant_Caspofungin'] == False)]
    V_Ani_Caspo_sequences = np.array2string(V_Ani_Caspo['Sequence'].values)

    V_Ani_Mica = df.loc[(df['Resistant_Anidulafungin'] == False) & (df['Resistant_Micafungin'] == False)]
    V_Ani_Mica_sequences = np.array2string(V_Ani_Mica['Sequence'].values)

    V_Caspo_Mica = df.loc[(df['Resistant_Caspofungin'] == False) & (df['Resistant_Micafungin'] == False)]
    V_Caspo_Mica_sequences = np.array2string(V_Caspo_Mica['Sequence'].values)

    if len(R_ALL_sequences) > 13:
        GenerateWebLogo(R_ALL_sequences, f'V3_{strain}_{locus}_all_R')

    if len(V_ALL_sequences) > 13:
        GenerateWebLogo(V_ALL_sequences, f'V3_{strain}_{locus}_all_V')

    if len(R_Ani_Caspo_sequences) > 13:
        GenerateWebLogo(R_Ani_Caspo_sequences, f'V3_{strain}_{locus}_Ani_Caspo_R')

    if len(R_Ani_Mica_sequences) > 13:
        GenerateWebLogo(R_Ani_Mica_sequences, f'V3_{strain}_{locus}_Ani_Mica_R')

    if len(R_Caspo_Mica_sequences) > 13:
        GenerateWebLogo(R_Caspo_Mica_sequences, f'V3_{strain}_{locus}_Caspo_Mica_R')

    if len(V_Ani_Caspo_sequences) > 13:
        GenerateWebLogo(V_Ani_Caspo_sequences, f'V3_{strain}_{locus}_Ani_Caspo_V')

    if len(V_Ani_Mica_sequences) > 13:
        GenerateWebLogo(V_Ani_Mica_sequences, f'V3_{strain}_{locus}_Ani_Mica_V')

    if len(V_Caspo_Mica_sequences) > 13:
        GenerateWebLogo(V_Caspo_Mica_sequences, f'V3_{strain}_{locus}_Caspo_Mica_V')

makeLogoPairPlots('BY4741', 'FKS1-HS1')