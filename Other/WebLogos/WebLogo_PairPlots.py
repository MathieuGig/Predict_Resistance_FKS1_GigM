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
    print(df)

    R_ALL = df.loc[(df['Resistant_Anidulafungin'] == True) & (df['Resistant_Caspofungin'] == True) & (df['Resistant_Micafungin'] == True)]
    R_ALL_values = R_ALL['Sequence'].values
    R_ALL_sequences = np.array2string(R_ALL_values)

    R_Ani_Caspo = df.loc[(df['Resistant_Anidulafungin'] == False) & (df['Resistant_Caspofungin'] == True) & (df['Resistant_Micafungin'] == True)]
    R_Ani_Caspo_values = R_Ani_Caspo['Sequence'].values
    R_Ani_Caspo_sequences = np.array2string(R_Ani_Caspo_values)

   # Que faire Que faire... tant de choix... tant de possibilites

    #GenerateWebLogo(R_ALL_sequences, f'{strain}_{locus}_all_resistant')

makeLogoPairPlots('BY4741', 'FKS1-HS1')