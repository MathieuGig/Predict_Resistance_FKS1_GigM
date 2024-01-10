import pandas as pd
import numpy as np
import re

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)


df1 = pd.read_csv('TableauFKS1.csv')
df2 = pd.read_csv('TableauFKS2.csv')

PathoDF1 = df1.loc[df1['Is_Human_Pathogen'] == True]

FKS1_Seq_from_pathogens = PathoDF1['Hotspot1'].values

TheSequences = np.array2string(FKS1_Seq_from_pathogens)

def GenerateWebLogo(sequences, name):
    x = sequences

    # Run some REGEX
    x = re.sub('[\[\]]', '', x)
    x = re.sub(' ', '', x)
    x = re.sub('\'\'', '\n', x)
    x = re.sub('[\'\[\]]', '', x)

    filename = f'WebLogoInput_{name}'

    f = open(filename, 'w')
    f.write(x)
    f.close()

    logo_name = f'logo_{name}'

    command = f'weblogo -f {filename} -o {logo_name} -F png -s large --errorbars NO'

    print(command)

#GenerateWebLogo(TheSequences, 'FKS1_FromPathogens')