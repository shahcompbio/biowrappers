'''
Created on Nov 9, 2015

@author: Andrew Roth
'''
from functools import wraps
from scipy.misc import logsumexp as log_sum_exp
from scipy.special import psi

import dask.array as da
import numpy as np
import scipy.stats as stats

def compute_e_log_dirichlet(x):
    return psi(x) - psi(np.expand_dims(x.sum(axis=-1), axis=-1))

def compute_e_log_q_dirichlet(x):
    return stats.dirichlet.entropy(x)

def compute_e_log_q_discrete(log_x):
    return np.sum(safe_multiply(np.exp(log_x), log_x))

def safe_multiply(x, y):
    return np.sign(x) * np.sign(y) * np.exp(np.log(np.abs(x)) + np.log(np.abs(y)))

def log_space_normalise(x):
    return x - log_sum_exp(x, axis=1)[:, np.newaxis]

log_sum_exp_chunk = da.chunk.keepdims_wrapper(log_sum_exp)

@wraps(log_sum_exp_chunk)
def log_sum_exp_da(a, axis=None, dtype=None, keepdims=False):
    return da.reductions.reduction(
        a, 
        log_sum_exp_chunk, 
        log_sum_exp_chunk, 
        axis=axis, 
        keepdims=keepdims,
        dtype='f8'
    )

def log_space_normalise_da(x):
    return x - log_sum_exp_da(x, axis=1)[:, np.newaxis]