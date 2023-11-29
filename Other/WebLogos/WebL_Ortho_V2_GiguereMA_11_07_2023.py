import pandas as pd
import re

df = pd.read_csv('epiortho.csv')
List = df['aaseqs'].values

row = 0

for i in List:

    strain = df.at[row, 'strain']
    locus = df.at[row, 'locus']
    compound = df.at[row, 'compound']
    resistance = df.at[row, 'resistance']
    filename = f'WebLogoInput_{strain}_{locus}_{compound}_{resistance}'

    x = re.sub('[\[\]]', '', i)
    x = re.sub(' ', '', x)
    x = re.sub('\'\'', '\n', x)
    x = re.sub('[\'\[\]]', '', x)

    f = open(filename, 'w')
    f.write(x)

    row += 1
    f.close()

    logo_name = f'logo_{strain}_{locus}_{compound}_{resistance}'
    command = f'weblogo -f {filename} -o {logo_name} -F png -s large --errorbars NO'

    print(command)