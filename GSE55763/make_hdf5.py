#!/usr/bin/env python3
import h5py
from tqdm import trange
import numpy as np

NUMROWS = 485577
NUMSAMPLES = 2711

with h5py.File('GSE55763.hdf5', 'w') as f:
    reads = f.create_dataset('reads', (NUMROWS, NUMSAMPLES, 2), dtype='int64')
    with open('GSE55763_unmethylated_methylated_signal_intensities.txt', 'r') as g:
        header = g.readline().split('\t')[1::3]
        header = [str(val.split(' ')[0]) for val in header]
        cpg_sites = []
        for i in trange(NUMROWS):
            line = g.readline()
            reads[i,:,0] = np.array(line.split('\t')[2::3]).astype(np.int64)
            reads[i,:,1] = np.array(line.split('\t')[1::3]).astype(np.int64)
            cpg_sites.append(line.split('\t')[0])
    sites = f.create_dataset('sites', data=header)
    smpls = f.create_dataset('smpls', data=cpg_sites)
