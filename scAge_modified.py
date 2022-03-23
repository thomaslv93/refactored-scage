import numpy as np
import pandas as pd

import scipy.stats as ss
from sklearn.linear_model import LinearRegression

import subprocess

from tqdm.contrib.concurrent import process_map

import multiprocessing
num_total_cores = multiprocessing.cpu_count()

# Takes a dataframe of training methylation data and returns sample names, cpg names, methylation levels, sample ages, and metadata
# Expects dataframe to have each row be a sample, and to have CpG site, age, and metadata in the columns
def _extract_data(train_df, metacols=["Tissue", "Strain", "Gender"]):
    # Extract basic information
    sample_names = list(train_df.index)                        
    cpg_names    = [x for x in train_df.columns if "chr" in x]
    metlev       = train_df[cpg_names].transpose().values
    ages         = train_df.loc[:, "Age"].values                             
    
    # Extract meta information
    metacols = list(train_df.columns.intersection(metacols))
    metacol_info = {}
    for col in metacols: metacol_info[col] = train_df[col].value_counts()
    
    return sample_names, cpg_names, metlev, ages, metacol_info


# Using methylation levels from a single CpG site from multiple samples and corresponding ages,
# perform linear regression from age to methylation levels
def _age_vs_metlev_regression(metlev, ages):
    # Deal only with non-nan methylation levels
    valid_metlevs = metlev[~np.isnan(metlev)]
    valid_ages = ages[~np.isnan(metlev)]

    # get Pearson correlation metrics                                           
    pr, pp = ss.pearsonr(valid_metlevs, valid_ages)                    

    # calculate linear regression (with traditional ordinary least squares (OLS))
    # Given an age, we will use this to estimate the methylation level.
    reg = LinearRegression().fit(valid_ages.reshape(-1, 1), valid_metlevs)           

    # isolate coefficient (slope) and intercept                                 
    coef = reg.coef_[0]                                                         
    intercept = reg.intercept_     

    # return output tuple                                                       
    ref_output = (pr, pp, coef, intercept)                   
    return ref_output    


# Run _age_vs_metlev_regression in parallel and format output
def _parallelize_regression(cpg_names, metlev, ages, n_cores=num_total_cores, chunksize=10_000):
    # Run _age_vs_metlev_regression in parallel
    results = process_map(_age_vs_metlev_regression,
                          metlev,
                          np.tile(ages, (metlev.shape[0], 1)),
                          max_workers = n_cores,
                          chunksize = chunksize) 
    
    # Final touches and formatting of output
    columns = ["PearsonR", "PearsonP", "Coef", "Intercept"]
    results = pd.DataFrame(results, columns=columns)
    results["ChrPos"] = cpg_names
    return results.set_index("ChrPos").dropna(axis=0)


# Extract data from the training dataframe, and run the parallelized regression algorithm
def construct_reference(train_df, n_cores=num_total_cores, chunksize=10_000):
    
    # Extract the data from the dataframe
    _, cpg_names, metlev, ages, _ = _extract_data(train_df)
    
    # Run the parallelized regression function on the data
    return _parallelize_regression(cpg_names, metlev, ages, n_cores=n_cores, chunksize=chunksize)


# Constants
PLB = .001 # probability lower bound (avoid zero probability)
PUB = .999 # Probability upper bound

# Helper functions

# Create column vector from numpy array
def _column_vectorify(vec): return np.expand_dims(vec.flatten(), axis=0).T

# Get a list of ages with regular intervals
def _age_gradient(min_age=-20, max_age=60, step=0.1): return np.arange(min_age, max_age + step, step).reshape(1,-1)

# Create a dictionary from a directory path
def _get_dict_from_dir(p): return {f.split(".tsv")[0]: os.path.join(p, f) for f in os.listdir(p) if f != ".ipynb_checkpoints"}

# Return the probability of agreement between two random events with probabilities p and q
def _probability_of_agreement(p, q): return p*q+(1-p)*(1-q)

def _probability_difference(p, q): return 1-np.abs(q-p)

# Project vector v onto the line defined by mx+b for x in [-inf,inf] and return magnitude
def _projection_magnitude(v, m, b):
    # dot products require one-dimensional vectors
    v, m, b = v.flatten(), m.flatten(), b.flatten()
    return np.dot((v-b), m)/np.dot(m, m)

# Get the log likelihood of overall methylation levels at each age
def _get_age_log_likelihoods(metlev,
                             slopes,
                             intercepts,
                             ages=_age_gradient(),
                             prob_func=_probability_of_agreement):
    # slopes and intercepts are column vectors; ages is a row vector
    # Use regression model to get the methylation estimate for each age step for each CpG site
    # Clip at probability lower bound (PLB) below and probability upper bound (PUB) above,
    # to preserve probabilistic interpretation and avoid
    # probabilities of 0 and 1...
    metlev, slopes, intercepts = _column_vectorify(metlev), _column_vectorify(slopes), _column_vectorify(intercepts)
    
    met_estimate = np.clip(slopes @ ages + intercepts, PLB, PUB)
    
    # Get the probability of similarity between actual methylation levels
    # and age-predicted methylation levels
    prob = prob_func(np.tile(metlev, (1, ages.shape[0])), met_estimate)
    # Return the sum of the log probabilities
    return np.sum(np.log(prob), axis=0)

# A version of _get_age_log_likelihoods which allows for weighting of probabilities
def _weighted_age_log_likelihoods(metlev,
                                  slopes,
                                  intercepts,
                                  weights,
                                  ages=_age_gradient(),
                                  prob_func=_probability_of_agreement):
    # slopes and intercepts are column vectors; ages is a row vector
    # Use regression model to get the methylation estimate for each age step for each CpG site
    # Clip at probability lower bound (PLB) below and probability upper bound (PUB) above,
    # to preserve probabilistic interpretation and avoid
    # probabilities of 0 and 1...
    metlev = _column_vectorify(metlev)
    slopes = _column_vectorify(slopes)
    intercepts = _column_vectorify(intercepts)
    weights = _column_vectorify(weights)
    
    met_estimate = np.clip(slopes @ ages + intercepts, PLB, PUB)
    
    # Get the probability of similarity between actual methylation levels
    # and age-predicted methylation levels
    prob = prob_func(np.tile(metlev, (1, ages.shape[0])), met_estimate)
    # Return the sum of the log probabilities
    return np.sum(weights * np.log(prob), axis=0)
    
# Get most likely age from log likelihoods
def _age_from_log_prob(metlev,
                       slopes,
                       intercepts,
                       parms=(_age_gradient(), _probability_of_agreement)):
    ages, prob_func = parms
    prob = _get_age_log_likelihoods(metlev, slopes, intercepts, ages, prob_func)
    return ages[:, np.argmax(prob)][0] # Select the age of maximum probability

# Get most likely age from log likelihoods
def _age_from_weighted_log_prob(metlev,
                                slopes,
                                intercepts,
                                weights=[],
                                parms=(_age_gradient(), _probability_of_agreement)):
    # Default is to weight all probabilities equally
    if len(weights)==0: weights = np.ones(metlev.shape[0])
    ages, prob_func = parms
    prob = _weighted_age_log_likelihoods(metlev, slopes, intercepts, weights, ages, prob_func)
    return ages[:, np.argmax(prob)][0] # Select the age of maximum probability

# Get confidence interval from log probabilities and ages
def _confidence_interval(log_probs, ages, ci=0.95):
    log_probs, ages = log_probs.flatten(), ages.flatten()
    
    assert log_probs.shape == ages.shape
    
    # Take exponential of log probabilities than normalize so that total probability is 1
    implicit_distribution = np.exp(log_probs)/np.sum(np.exp(log_probs))
    
    # Location of most probable age
    peak = np.argmax(implicit_distribution)
    
    # Get cumulative values out from the point of largest probability
    leftside = np.cumsum(np.flip(implicit_distribution[0:peak+1])) - implicit_distribution[peak]/2
    rightside = np.cumsum(implicit_distribution[peak:-1]) - implicit_distribution[peak]/2
    
    # Get offset for half of confidence interval leftward and rightward
    left_offset = np.argwhere(leftside < ci/2)[-1][0]
    right_offset = np.argwhere(rightside < ci/2)[-1][0]
    
    return ages[peak - left_offset], ages[peak + right_offset]

# Distance-based alternative to probabilistic approach; similar results
def _age_from_projection(metlev,
                         slopes,
                         intercepts,
                         parms=()):
    return _projection_magnitude(metlev, slopes, intercepts)

# Take df input with a 'PearsonR' column and selects
# return view not copy
def _quantile_select(df, parm=0.99):
    quantile = df['PearsonR'].abs().quantile(q=parm)
    idx = df[df['PearsonR'].abs() >= quantile].index
    return df.loc[idx]

def _number_select(df, parm=500): return df.loc[df['PearsonR'].abs().nlargest(parm).index]
    
def _cutoff_select(df, parm=0.7): return df.loc[df['PearsonR'].abs() >= parm]

# Extract CpG names, methylation levels, slopes, and intercepts from input data and the reference dataframe
def _filter_cpgs(in_df, ref_df, selection_function=_quantile_select, parm=0.99):
    # Keep only the CpG sites which the data has in common with the reference
    in_df = pd.concat([ref_df, in_df], axis=1, join='inner')
    
    # Keep only the top CpG sites sorted by Pearson's R
    in_df = selection_function(in_df, parm)
    
    # Extract and return the important data
    cpg_names  = in_df.index.values
    metlev     = in_df['MetLev'].values.flatten().astype(bool).reshape(-1,1)
    slopes     = in_df['Coef'].values.flatten().reshape(-1,1)
    intercepts = in_df['Intercept'].values.flatten().reshape(-1,1)
    
    return cpg_names, metlev, slopes, intercepts

def _get_age(in_df,
             ref_df,
             age_func=_age_from_log_prob,
             age_parms=(_age_gradient(), _probability_of_agreement),
             sel_func=_quantile_select,
             sel_parm=0.99):
    _, metlev, slopes, intercepts = _filter_cpgs(in_df, ref_df, sel_func, sel_parm)
    return age_func(metlev, slopes, intercepts, age_parms)

# Run the _likeliest_age function in parallel
def run_scAge(in_dfs,    # dictionary of input dataframes ({samplename: dateframe, ...})
              ref_df,
              age_func=_age_from_log_prob,
              age_parms=(_age_gradient(), _probability_of_agreement),
              sel_func=_quantile_select,
              sel_parm=0.99,
              n_cores = num_total_cores,
              chunksize = 1):
    
    # Run _likeliest_age in parallel
    results = process_map(_get_age,
                          list(in_dfs.values()),   # single sample methylation dataframe
                          len(in_dfs)*[ref_df],    # reference dataframe
                          len(in_dfs)*[age_func],  # age function
                          len(in_dfs)*[age_parms], # parameters for age function
                          len(in_dfs)*[sel_func],  # cpg selection function
                          len(in_dfs)*[sel_parm],  # parameters for cpg selection function
                          max_workers = n_cores,
                          chunksize = chunksize)
    
    return pd.DataFrame(results, columns=['PredictedAge'], index=in_dfs.keys())