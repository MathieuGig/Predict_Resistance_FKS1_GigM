import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import glob
import os

def Epistasis(condition):
    # Concatenate all dataframes
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

    ShortData = master[['Nmut_codons', 'nt_seq', 'alt_codons', 'aa_pos', 'median_s', 'strain', 'locus', 'compound']]
    ShortData = ShortData[ShortData['strain'] == condition]

    SingleMutations = ShortData[ShortData['Nmut_codons'] == 1.0]

    MultipleMutations = ShortData[ShortData['Nmut_codons'] > 1.0]
    MultipleMutations = MultipleMutations.sort_values(by='Nmut_codons')
    MultipleMutations = MultipleMutations.sort_values(by='aa_pos')

    AdditivityScore = list()
    OrthologScore = list()
    ListLocus = []
    ListCompound = []

    ListeSequencesOrtholog = list(set(MultipleMutations['nt_seq']))

    AAtoCodon = {'S': ['tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'], 'L': ['tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'],
                 'C': ['tgt', 'tgc'], 'W': ['tgg'], 'E': ['gaa', 'gag'], 'D': ['gat', 'gac'],
                 'P': ['cct', 'ccc', 'cca', 'ccg'],
                 'V': ['gtt', 'gtc', 'gta', 'gtg'], 'N': ['aat', 'aac'], 'M': ['atg'], 'K': ['aaa', 'aag'],
                 'Y': ['tat', 'tac'],
                 'I': ['att', 'atc', 'ata'], 'Q': ['caa', 'cag'], 'F': ['ttt', 'ttc'],
                 'R': ['cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'], 'T': ['act', 'acc', 'aca', 'acg'],
                 '*': ['taa', 'tag', 'tga'],
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
        ListLocus.append(Mutations['locus'].values[0])
        ListCompound.append(Mutations['compound'].values[0])

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
                print('no synonymous codon found')
                break

            selCo = selCo['median_s'].values[0]

            additivity.append(selCo)

        additivity = sum(additivity)
        AdditivityScore.append(additivity)

    # Now do graphs
    x = AdditivityScore
    y = OrthologScore
    z = ListLocus
    w = ListCompound
    graphdf = pd.DataFrame(list(zip(x, y, z, w)), columns=['Additivity', 'Ortholog', 'Locus', 'Compound'])
    #print(graphdf)
    sns.set_theme(style='whitegrid')

    compounds = ['caspofungin', 'micafungin', 'anidulafungin', 'none']
    SortLocus = ['FKS1-HS1', 'FKS1-HS2', 'FKS2-HS1', 'FKS2-HS2']

    graphdf = graphdf[(graphdf.Compound.isin(compounds))]

    fig = sns.relplot(graphdf,
                      x='Additivity', y='Ortholog',
                      row='Locus', col='Compound', col_order=compounds, row_order=SortLocus,
                      height=2.5, aspect=1,
                      )

    #fig.set(xlim=(-5, 10), ylim=(-5, 5))

    fig.set_titles(row_template='{row_name}', col_template='{col_name}')
    fig.set(xlabel='Additivity', ylabel='Ortholog Sel_coef')

    fig.fig.subplots_adjust(top=.88)
    fig.fig.suptitle(condition)

    plt.show()



if __name__ == '__main__':
    Epistasis('BY4741')
    Epistasis('R1158')

