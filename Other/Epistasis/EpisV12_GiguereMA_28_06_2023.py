import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import glob
import os
import numpy as np

# Voir pandas.tuple


pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

def epistasis():
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

    per_mut = master[master.seq_type.isin(['single', 'ortho'])].pivot_table(
        index=['strain', 'locus', 'compound', 'nt_seq', 'seq_type', 'median_s'],
        columns='mutated_codon', values=['aa_pos', 'alt_codons'],
        aggfunc='first').reset_index()


    pd.set_option('mode.chained_assignment', None)

    #Orthologs
    per_mut_ortholog = per_mut[per_mut.seq_type == 'ortho']
    #per_mut_ortholog = per_mut_ortholog.head()
    per_mut_ortholog.columns = [x[0] for x in per_mut_ortholog.columns[:-18]] + [f"{x[0]}_{int(x[1])}" for x in per_mut_ortholog.columns[-18:]]
    per_mut_ortholog.rename(columns={'median_s_1':'s_double'}, inplace=True)
    #per_mut_ortholog.drop(columns=['seq_type', 'median_s_2'], inplace=True)

    #per_mut_ortholog.columns = ['strain', 'locus', 'compound', 'nt_seq', 'seq_type', 'median_s',
                                #'aa_pos_1', 'aa_pos_2', 'aa_pos_3',  'aa_pos_4', 'aa_pos_5', 'aa_pos_6', 'aa_pos_7', 'aa_pos_8', 'aa_pos_9',
                                #'alt_codons_1', 'alt_codons_2', 'alt_codons_3', 'alt_codons_4', 'alt_codons_5', ]
    print(per_mut_ortholog)
    #print(per_mut_ortholog[('aa_pos', 7.0)])

    per_mut_single = per_mut[per_mut.seq_type == 'single']
    per_mut_single.rename(columns={'median_s': 's_single'}, inplace=True)
    per_mut_single.columns = [x[0] for x in per_mut_single.columns[:-18]] + [f"{x[0]}_{int(x[1])}" for x in per_mut_single.columns[-18:]]
    #per_mut_single = per_mut_single.dropna(axis='columns')

    print(per_mut_single)

    per_mut_concat = pd.merge(left=per_mut_ortholog, right=per_mut_single, how='outer', indicator='location_1', suffixes=(None, '1'),
                              on=['strain', 'locus', 'compound', 'aa_pos_1', 'alt_codons_1'])

    print(per_mut_concat)

    for i in [2, 3, 4, 5, 6, 7, 8, 9]:
        print(f'location_{i}')
        per_mut_concat_merged = pd.merge(left=per_mut_concat, right=per_mut_single, how='outer', indicator=f'location_{i}', suffixes=(None, f'_{i}'),
                                  on=['strain', 'locus', 'compound', f'aa_pos_{i}', f'alt_codons_{i}'])

    #print(per_mut_concat.groupby('location_9').size())
    print(per_mut_concat_merged)



    #Graphs
    # See V10



if __name__ == '__main__':
    epistasis()