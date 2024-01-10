## Title: GenerateWebLogo
## Author: Mathieu Giguere
## Brief: Uses regex to format a list of strings into a file with 1 string per line,
##        then prints a command for "weblogo3" that can be copy-pasted into a terminal.

import pandas as pd
import numpy as np
import re

pd.set_option('display.max_columns', None)
pd.set_option('display.expand_frame_repr', False)
pd.set_option('max_colwidth', None)

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