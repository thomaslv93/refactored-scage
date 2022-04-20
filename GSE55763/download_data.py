#!/usr/bin/env python3
import os
import numpy as np
from tqdm import trange

# Download the appropriate dataset
ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55763/suppl/GSE55763_unmethylated_methylated_signal_intensities.txt.gz'
filename = 'GSE55763_unmethylated_methylated_signal_intensities.txt.gz'
os.system(f'wget {ftp_url}')
os.system('gunzip *gz')

NUM_ROWS = 485_578

# Convert data to binary files where rows and columns can be quickly accessed
with open('GSE55763_unmethylated_methylated_signal_intensities.txt', 'r') as f:
    with open('unmethylated_bytes', 'wb') as ub:
        with open('methylated_bytes', 'wb') as mb:
            line = f.readline() # Pass over header
            for i in trange(NUM_ROWS):
                if not line: break
                line = f.readline()
                unm = np.array(line.split('\t')[1::3]).astype(int)
                met = np.array(line.split('\t')[2::3]).astype(int)
                for val in unm:
                    ub.write(int(val).to_bytes(4, 'little'))
                for val in met:
                    mb.write(int(val).to_bytes(4, 'little'))
                i += 1
