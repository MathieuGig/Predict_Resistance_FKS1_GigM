import pandas as pd
import re

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)


df = pd.read_csv('epiortho.csv')
List = df['aaseqs'].values
strain = df['strain'].values
print(strain)
for i in List:
    #print(df[df['aaseqs'] == i].values)
    x = re.sub(' ', '', i)
    x = re.sub('\'\'', '\n', x)
    x = re.sub('[\'\[\]]', '', x)
    #print(x)

print('the command')