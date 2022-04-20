#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os
from tqdm.contrib.concurrent import process_map


import sys
sys.path.append('..')
import scAge_modified as scAge

FOLDERPATH = os.path.dirname(__file__) + '/'

CPG_SITENAMES = list(pd.read_csv(FOLDERPATH+'signal_sitenames.csv').values.flatten())
SAMPLE_NAMES = list(pd.read_csv(FOLDERPATH+'signal_samplenames.csv').values.flatten())
BYTE_SIZE = 4
NUM_SAMPLES = len(SAMPLE_NAMES)
NUM_CPG_SITES = len(CPG_SITENAMES)

all_meta = pd.read_csv(FOLDERPATH+'all-meta.csv', index_col=0)
AGES = (all_meta.loc[SAMPLE_NAMES][['Age']]).values.flatten()

CHUNK_SIZE = NUM_CPG_SITES//480 + 1
CHUNKS = [list(range(i*CHUNK_SIZE, min((i+1)*CHUNK_SIZE, NUM_CPG_SITES))) for i in range((NUM_CPG_SITES//CHUNK_SIZE)+1)]

def _generate_reference(row_idxs, col_idxs):
    with open(FOLDERPATH+'methylated_bytes', 'rb') as met:
        with open(FOLDERPATH+'unmethylated_bytes', 'rb') as unm:
            parameters = []
            for row_idx in row_idxs:
                met_vals = []
                unm_vals = []
                
                # Build row
                for col in col_idxs:
                    met.seek((NUM_SAMPLES*row_idx+col)*BYTE_SIZE,0)
                    unm.seek((NUM_SAMPLES*row_idx+col)*BYTE_SIZE,0)
                    met_val = met.read(BYTE_SIZE)
                    unm_val = unm.read(BYTE_SIZE)
            
                    try: met_vals.append(np.frombuffer(met_val, dtype='int32')[0])
                    except: continue
                    try: unm_vals.append(np.frombuffer(unm_val, dtype='int32')[0])
                    except: continue
                    
                met_vals = np.array(met_vals)
                unm_vals = np.array(unm_vals)
                with np.errstate(divide='ignore', invalid='ignore'):
                    metlev = met_vals / (met_vals + unm_vals)
                parameters.append(scAge._age_vs_metlev_regression(metlev, AGES[col_idxs]))
    return parameters

def generate_reference():
    parms = process_map(_generate_reference,
	                    CHUNKS,
	                    [range(70)]*len(CHUNKS),
	                    max_workers=60,
	                    chunksize=1)
    ewas_ref_df = []
    for parm_chunk in parms:
        ewas_ref_df.extend(parm_chunk)
    ewas_ref_df = pd.DataFrame(ewas_ref_df, columns=['PearsonR', 'PearsonP', 'Coef', 'Intercept'], index=CPG_SITENAMES)
    ewas_ref_df.to_csv('./ewas_ref_df.csv')
    return ewas_ref_df 

def get_ages():
    return AGES

if __name__ == "__main__": generate_reference()
