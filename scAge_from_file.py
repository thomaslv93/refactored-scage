#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm
import scAge_modified as scAge
import heapq

# Convert the reference dataframe to a dictionary for rapid access
# Believe it or not, this is significantly faster than a pandas dataframe
# for simple data retrieval...
def ref_to_dict(ref_df):
    return dict(zip(ref_df.index, 
                    list(zip(ref_df['PearsonR'],
                             ref_df['PearsonP'],
                             ref_df['Coef'],
                             ref_df['Intercept']))))


# Get attributes from a reference dictionary.
# Return an empty list if the site is not in the dictionary.
def get_attribs(name, ref_dict): 
    try: attribs = ref_dict[name]
    except: attribs = []
    return tuple(attribs)


# A class representing a CpG site.
class CpGSite:
    def __init__(self, name, metlev, met, unm, pearsonr, pearsonp, coef, intercept):
        self.name = name      # CpG site name
        self.metlev = metlev  # Methylation level (in [0,1])
        self.met = met        # Number of methylated reads
        self.unm = unm        # Number of unmethylated reads
        
        # Pearson r correspondence of site with age according to reference data
        self.pearsonr = pearsonr
        self.pearsonp = pearsonp # p-value of Pearson's r
        
        # Coefficient of linear regression with age according to reference data
        self.coef = coef
        # Intercept of linear regression with age according to reference data
        self.intercept = intercept
    
    # Less than method, comparing absolute Pearson r
    def __lt__(self, other):
        return abs(self.pearsonr) < abs(other.pearsonr)
    
    # For easy comparison, print name and Pearson r value
    def __repr__(self):
        return f'({self.name}, {self.pearsonr})'


def naive_sampling_update(val, idx, readnum, selected):
    new_val = 0 # initialize the new count
    rem = val # track the remainder
    while True:
        if idx >= len(selected): break
        diff = selected[idx] - readnum
        if rem >= diff: # Increment new_met when we cross sample
            new_val += 1 
            rem = rem - diff # Update remainder
            readnum += diff # Increment read number
            idx += 1 # Increment sample idx
        else:
            readnum += rem
            break
    return new_val, idx, readnum


def count_reads(cpg_parms_getter):
    cumsum = 0
    for parms in cpg_parms_getter.yielder():
        cumsum += np.sum(parms[1:-1]).astype(int)
    return int(cumsum)


# This functions takes a context manager/iterator, a reference dictionary,
# a number of subsamples, and CpG parameter and returns the most age-correlated
# CpG sites according to these parameters
def read_initial_sites(cpg_parms_getter,
                       ref_dict,
                       nsubsamples=100_000,
                       ncpgs=10_000,
                       seed=1):
    
    rng = np.random.default_rng(seed)
    total_reads = count_reads(cpg_parms_getter)
    selected = np.sort(rng.choice(total_reads, nsubsamples, replace=False))

    heap = []
    readnum = 0
    idx = 0
    for cpg_parms in tqdm(cpg_parms_getter):
        name, met, unm, _ = cpg_parms # Throw away age value
        
        # Subsample to get a new count for unmethylated and methylated sites
        new_met, idx, readnum = naive_sampling_update(met, idx, readnum, selected)
        new_unm, idx, readnum = naive_sampling_update(unm, idx, readnum, selected)

        new_cov = new_met + new_unm
        
        if new_cov != 0:
            new_metlev = new_met / new_cov

            # Create a site object which can be compared to other
            # site objects
            attribs = get_attribs(name, ref_dict)
            if len(attribs)==0: 
                continue # if not in reference, skip

            site = CpGSite(name, new_metlev, new_met, new_unm, *attribs)

            # If the heap has reached its desired length...
            if len(heap) == ncpgs:
                # Only add site if bigger than current lowest site
                if heap[0] < site:
                    heapq.heappop(heap) # Remove lowest site
                    heapq.heappush(heap, site) # Add new site
            # If heap is less than desired length, add site...
            else:
                heapq.heappush(heap, site)

    # Get sites in ascending order of correlation
    sites = []
    while True:
        if len(heap)==0:
            break
        sites.append(heapq.heappop(heap))
        
    return sites

def age_from_sites(sites):
    # Perform scAge on selected CpG sites
    metlevs = np.array([site.metlev for site in sites]).reshape(-1,1)
    slopes = np.array([site.coef for site in sites]).reshape(-1,1)
    intercepts = np.array([site.intercept for site in sites]).reshape(-1,1)

    return scAge._age_from_log_prob(metlevs,
                                    slopes,
                                    intercepts,
                                    parms=(scAge._age_gradient(), scAge._probability_of_agreement))
