import pandas as pd

def print_hi(name):
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    a = pd.read_csv('DMS-main/processed_data/BY4741_FKS1-HS1_single_ortho_anidulafungin/selcoeff_all_libraries.csv')
    b = a[['Nmut_codons', 'aa_seq', 'alt_aa', 'median_s']]
    print(b['Nmut_codons'] == 1.0)
    c = b[b['Nmut_codons'] == 1.0]
    print(c)
    # b = a.groupby(['median_s', 'aa_seq', 'alt_aa'])
    # print(b)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
