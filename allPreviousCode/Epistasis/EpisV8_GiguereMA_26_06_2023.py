import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import glob
import os
from statistics import median

def Epistasis(condition, gene, hotspot, drug):
    FileName = 'DMS-main/processed_data/' + condition + '_' + gene + '-' + hotspot + ''

if __name__ == '__main__':

    #Concatenate all dataframes
    list_df = []
    for d in glob.glob('DMS-main/processed_data/*'):
        cset = os.path.basename(d).split('_')
        cset_name = '_'.join(cset)
        strain, locus, pool_type, compound = [cset[i] for i in (0, 1, 2, -1)]

        df = pd.read_csv(f'{d}/selcoeff_all_libraries.csv', index_col=0)
        df['strain'] = strain
        df['locus'] = locus
        df['pool_type'] = pool_type
        df['compound'] = compound

        if pool_type != 'double':
            list_df.append(df)

    master = pd.concat(list_df, ignore_index=True)
    master.groupby(['strain', 'locus',  # Per input
                    ])[['nt_seq']].nunique().reset_index()

#######################################################################

    TheData = pd.read_csv('DMS-main/processed_data/BY4741_FKS1-HS1_single_ortho_anidulafungin/selcoeff_all_libraries.csv')
    #TheData = master

    ShortData = TheData[['Nmut_codons', 'nt_seq', 'alt_codons', 'aa_pos', 'median_s']]
    #ShortData = TheData[['Nmut_codons', 'nt_seq', 'alt_codons', 'aa_pos', 'median_s', 'strain', 'locus', 'compound']]

    SingleMutations = ShortData[ShortData['Nmut_codons'] == 1.0]

    MultipleMutations = ShortData[ShortData['Nmut_codons'] > 1.0]
    MultipleMutations = MultipleMutations.sort_values(by='Nmut_codons')
    MultipleMutations = MultipleMutations.sort_values(by='aa_pos')

    AdditivityScore = list()
    OrthologScore = list()

    ListeSequencesOrtholog = list(set(MultipleMutations['nt_seq']))

    AAtoCodon = {'S': ['tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'], 'L': ['tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'],
     'C': ['tgt', 'tgc'], 'W': ['tgg'], 'E': ['gaa', 'gag'], 'D': ['gat', 'gac'], 'P': ['cct', 'ccc', 'cca', 'ccg'],
     'V': ['gtt', 'gtc', 'gta', 'gtg'], 'N': ['aat', 'aac'], 'M': ['atg'], 'K': ['aaa', 'aag'], 'Y': ['tat', 'tac'],
     'I': ['att', 'atc', 'ata'], 'Q': ['caa', 'cag'], 'F': ['ttt', 'ttc'],
     'R': ['cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'], 'T': ['act', 'acc', 'aca', 'acg'], '*': ['taa', 'tag', 'tga'],
     'A': ['gct', 'gcc', 'gca', 'gcg'], 'G': ['ggt', 'ggc', 'gga', 'ggg'], 'H': ['cat', 'cac']}

    CodonToAA = {
        'tca': 'S', 'tcc': 'S', 'tcg': 'S', 'tct': 'S',
        'ttc': 'F', 'ttt': 'F',
        'tta': 'L', 'ttg': 'L',
        'tac': 'Y', 'tat': 'Y',
        'taa': '*', 'tag': '*',
        'tgc': 'C', 'tgt': 'C',
        'tga': '*',
        'tgg': 'W',
        'cta': 'L', 'ctc': 'L', 'ctg': 'L', 'ctt': 'L',
        'cca': 'P', 'ccc': 'P', 'ccg': 'P', 'cct': 'P',
        'cac': 'H', 'cat': 'H',
        'caa': 'Q', 'cag': 'Q',
        'cga': 'R', 'cgc': 'R', 'cgg': 'R', 'cgt': 'R',
        'ata': 'I', 'atc': 'I', 'att': 'I',
        'atg': 'M',
        'aca': 'T', 'acc': 'T', 'acg': 'T', 'act': 'T',
        'aac': 'N', 'aat': 'N',
        'aaa': 'K', 'aag': 'K',
        'agc': 'S', 'agt': 'S',
        'aga': 'R', 'agg': 'R',
        'gta': 'V', 'gtc': 'V', 'gtg': 'V', 'gtt': 'V',
        'gca': 'A', 'gcc': 'A', 'gcg': 'A', 'gct': 'A',
        'gac': 'D', 'gat': 'D', 'gaa': 'E', 'gag': 'E',
        'gga': 'G', 'ggc': 'G', 'ggg': 'G', 'ggt': 'G'
    }

    for sequence in ListeSequencesOrtholog:

        Mutations = MultipleMutations.loc[MultipleMutations['nt_seq'] == sequence]

        OrthologScore.append(Mutations['median_s'].values[0])

        pos = list(Mutations.loc[:, 'aa_pos'])
        alt = list(Mutations.loc[:, 'alt_codons'])

        additivity = list()

        for i in range(len(pos)):
            selCo = SingleMutations.loc[(SingleMutations['aa_pos'] == pos[i]) & (SingleMutations['alt_codons'] == alt[i])]

            if selCo.empty:
                #Remplacer par codons synonymes.
                AA = CodonToAA.get(alt[i])
                Synonymous = AAtoCodon.get(AA)
                #selSyno = list()

                for j in range(len(Synonymous)):
                    selCo = SingleMutations.loc[(SingleMutations['aa_pos'] == pos[i]) & (SingleMutations['alt_codons'] == Synonymous[j])]
                    #print(selCo.get('median_s'))

                    if len(selCo) != 0:
                        break

            if selCo.empty:
                break

            selCo = selCo['median_s'].values[0]

            additivity.append(selCo)

        additivity = sum(additivity)
        AdditivityScore.append(additivity)

    # Now do graphs
    x = AdditivityScore
    y = OrthologScore
    df = pd.DataFrame(list(zip(x, y)), columns=['Additivity', 'Ortholog'])
    #print(df)
    sns.set_theme(style='white')
    fig = sns.relplot(data=df, x="Additivity", y="Ortholog")
    plt.show()