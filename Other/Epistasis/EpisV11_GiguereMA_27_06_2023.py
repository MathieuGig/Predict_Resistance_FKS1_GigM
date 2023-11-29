import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import glob
import os
import numpy as np

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

    per_mut = master[master.seq_type.isin(['single', 'ortho'])].pivot_table(
        index=['strain', 'locus', 'compound', 'nt_seq', 'seq_type'],
        columns='mutated_codon', values=['aa_pos', 'alt_codons', 'median_s'],
        aggfunc='first').reset_index()


    pd.set_option('mode.chained_assignment', None)

    #Orthologs
    per_mut_ortholog = per_mut[per_mut.seq_type == 'ortho']
    print(per_mut_ortholog)

    per_mut_single = per_mut[per_mut.seq_type == 'single']
    per_mut_single.rename(columns={'median_s': 's_single'}, inplace=True)
    per_mut_single.columns = [x[0] for x in per_mut_single.columns[:-6]] + [f"{x[0]}_{int(x[1])}" for x in
                                                                            per_mut_single.columns[-6:]]
    per_mut_single = per_mut_single.dropna(axis='columns')

    print(per_mut_single)

    #Graphs



if __name__ == '__main__':
    Epistasis('BY4741')
    #Epistasis('R1158')