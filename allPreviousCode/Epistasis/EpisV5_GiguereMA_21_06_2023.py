import pandas as pd
from statistics import median

def Epistasis(condition, gene, hotspot, drug):
    FileName = 'DMS-main/processed_data/' + condition + '_' + gene + '-' + hotspot + ''

if __name__ == '__main__':
    TheData = pd.read_csv('DMS-main/processed_data/BY4741_FKS1-HS1_single_ortho_anidulafungin/selcoeff_all_libraries.csv')

    ShortData = TheData[['Nmut_codons', 'nt_seq', 'alt_codons', 'aa_pos', 'median_s']]

    SingleMutations = ShortData[ShortData['Nmut_codons'] == 1.0]

    MultipleMutations = ShortData[ShortData['Nmut_codons'] > 1.0]
    MultipleMutations = MultipleMutations.sort_values(by='Nmut_codons')
    MultipleMutations = MultipleMutations.sort_values(by='aa_pos')

    ListeSequencesOrtholog = list(set(MultipleMutations['nt_seq']))

    for sequence in ListeSequencesOrtholog:
        print(sequence)
        Mutations = MultipleMutations.loc[MultipleMutations['nt_seq'] == sequence]
        print(Mutations)
        pos = list(Mutations.loc[:, 'aa_pos'])
        alt = list(Mutations.loc[:, 'alt_codons'])
        print(pos)
        print(alt)
        additivity = list()
        for i in range(len(pos)):
            selCo = SingleMutations.loc[(SingleMutations['aa_pos'] == pos[i]) & (SingleMutations['alt_codons'] == alt[i])]
            selCo = selCo['median_s'].values[0]
            #selCo = list(selCo.loc[:, 'median_s'])
            #selCo = selCo.get(['median_s'])
            print(selCo)
            additivity.append(selCo)
            #print(additivity)
        additivity = sum(additivity)
        #print(additivity)

        #df.loc[df['column_name'] == some_value]