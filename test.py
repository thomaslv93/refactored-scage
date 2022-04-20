#!/usr/bin/env python3
from scAge_from_file import read_initial_sites, ref_to_dict, age_from_sites
import pandas as pd

import sys
sys.path.append('./GSE55763')
from generate_reference import generate_reference, get_ages
from iterator import ewas_cm_it

ewas_ref_df = pd.read_csv('./GSE55763/ewas_ref_df.csv', index_col=0)
AGES = get_ages()

seed = 1
nsubsamples = 100_000
ncpgs = 10_000
ref_dict = ref_to_dict(ewas_ref_df)

if len(sys.argv) > 1:
    try:
        arg = [int(sys.argv[1])]
    
    except ValueError:
        print('arg must be an int')
else:
    arg = [0]

with ewas_cm_it(arg) as e:
    sites = read_initial_sites(cpg_parms_getter=e,
                               ref_dict=ref_dict,
                               nsubsamples=nsubsamples,
                               ncpgs=ncpgs,
                               seed=seed)

age = age_from_sites(sites)
print(f'predicted age: {age}')
print(f'actual age: {AGES[0]}')
