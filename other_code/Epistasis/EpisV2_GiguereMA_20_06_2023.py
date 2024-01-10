import pandas as pd

def Epistasis(condition, gene, hotspot, drug):
    FileName = 'DMS-main/processed_data/' + condition + '_' + gene + '-' + hotspot + ''

if __name__ == '__main__':
    TheData = pd.read_csv('DMS-main/processed_data/BY4741_FKS1-HS1_single_ortho_anidulafungin/selcoeff_all_libraries.csv')

    ShortData = TheData[['Nmut_codons', 'aa_seq', 'alt_aa', 'aa_pos', 'median_s']]

    SingleMutations = ShortData[ShortData['Nmut_codons'] == 1.0]
    #print(SingleMutations)

    MultipleMutations = ShortData[ShortData['Nmut_codons'] > 1.0]
    MultipleMutations = MultipleMutations.sort_values(by='Nmut_codons')
    MultipleMutations = MultipleMutations.sort_values(by='aa_pos')
    #print(MultipleMutations)

    ListeSequencesOrtholog = list(set(MultipleMutations['aa_seq']))

    for sequence in ListeSequencesOrtholog:
        #print(sequence)
        Mutations = MultipleMutations.loc[MultipleMutations['aa_seq'] == sequence]
        #print(Mutations)
        pos = Mutations.loc[:, 'aa_pos']
        alt = Mutations.loc[:, 'alt_aa']
        #miam = SingleMutations.loc[(SingleMutations['aa_pos'] == pos) & (SingleMutations['alt_aa'] == alt)]
        #apply sum to miam


        #df.loc[df['column_name'] == some_value]