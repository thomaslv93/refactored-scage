#! /usr/bin/env python3
import numpy as np
import pandas as pd
import os

FOLDERPATH = os.path.dirname(__file__) + '/'

def sliceify(func, slce):
    start, stop, step = slce.start, slce.stop, slce.step
    if not slce.step: step = 1
    return [func(step * i + start) for i in range((stop - start) // step)]

# Define a context manager and iterator for the Epigenome-wide association study
# (EWAS) in dataset GSE55763
class ewas_cm_it():
    def __init__(self, col_idxs=[]):
        self.CPG_SITES = pd.read_csv(FOLDERPATH+'signal_sitenames.csv').values.flatten()
        self.SAMPLES = pd.read_csv(FOLDERPATH+'signal_samplenames.csv').values.flatten()
        self.AGES = pd.read_csv(FOLDERPATH+'all-meta.csv').set_index('Sample_Name')[['Age']].loc[self.SAMPLES].values.flatten()
        self.NUM_CPG_SITES = self.CPG_SITES.shape[0]
        self.NUM_SAMPLES = self.SAMPLES.shape[0]
        self.BYTE_SIZE = 4
        
        self.index = 0
        
        self._met_file_name = FOLDERPATH+'methylated_bytes'
        self._met_file_mode = 'rb'
        self._unm_file_name = FOLDERPATH+'unmethylated_bytes'
        self._unm_file_mode = 'rb'

        if not col_idxs: self.col_idxs = np.array(range(self.NUM_SAMPLES))
        else: self.col_idxs = np.array(col_idxs)

        if not np.issubdtype(self.col_idxs.dtype, np.integer):
            self.col_idxs = self._convert_samples_to_indices(col_idxs)
        
    # helper function to get individual index
    def _get_vals_helper(self, key, file_object):
        if key < 0: key = key + len(self) # Deal with negative indices
        if key >= len(self) or key < 0: raise IndexError # Deal with out of bounds indices
        vals = []
        for col in self.col_idxs:
            file_object.seek((self.NUM_SAMPLES*key+col)*self.BYTE_SIZE,0)
            result = file_object.read(self.BYTE_SIZE)
            if not result: raise IndexError
            try: vals.append(np.frombuffer(result, dtype='int32')[0])
            except: vals.append(np.array(0, dtype=np.int32))
        return np.array(vals)
    
    def _convert_cpg_to_index(self, cpg):
        return np.where(self.CPG_SITES==cpg)[0][0]
    
    def _convert_samples_to_indices(self, sample_list):
        return [np.where(self.SAMPLES==smpl)[0][0] for smpl in sample_list]
    
    def _get_vals(self, key):
        name = self.CPG_SITES[key]
        met_val = self._get_vals_helper(key, self._met_file)
        unm_val = self._get_vals_helper(key, self._unm_file)
        ages = self.AGES[self.col_idxs]
        return [name, met_val, unm_val, ages]
        
    # Indexing and slicing
    def __getitem__(self, key):
        if isinstance(key, slice):
            return np.array(sliceify(self._get_vals, key))
        try: key = int(key)
        except: key = self._convert_cpg_to_index(key)
        return self._get_vals(key)
    
    # Faster way to iterate
    def yielder(self):
        # Initialize files
        self._met_file.seek(0, 0)
        self._unm_file.seek(0, 0)
        for row in range(self.NUM_CPG_SITES):
            name = self.CPG_SITES[row]
            met_vals = []
            unm_vals = []
            for col in self.col_idxs:
                self._met_file.seek((self.NUM_SAMPLES*row+col)*self.BYTE_SIZE,0)
                self._unm_file.seek((self.NUM_SAMPLES*row+col)*self.BYTE_SIZE,0)
                result_met, result_unm = self._met_file.read(self.BYTE_SIZE), self._unm_file.read(self.BYTE_SIZE)
                if not result_met: raise StopIteration
                if not result_unm: raise StopIteration
                try: met_vals.append(np.frombuffer(result_met, dtype='int32')[0])
                except: met_vals.append(np.array(0, dtype=np.int32))
                try: unm_vals.append(np.frombuffer(result_unm, dtype='int32')[0])
                except: unm_vals.append(np.array(0, dtype=np.int32))
            
            yield [name, np.array(met_vals), np.array(unm_vals)]    
        
    # Make an iterator
    def __iter__(self):
        self.index = 0
        return self
    
    def __next__(self):
        if self.index < len(self):
            ret = self[self.index]
            self.index += 1
            return ret
        else:
            raise StopIteration
            
    # Length
    def __len__(self):
        return self.NUM_CPG_SITES
            
    # Context management
    def __enter__(self):
        self._met_file = open(self._met_file_name, self._met_file_mode)
        self._unm_file = open(self._unm_file_name, self._unm_file_mode)
        return self
    
    # Context management
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self._met_file.close()
        self._unm_file.close()
