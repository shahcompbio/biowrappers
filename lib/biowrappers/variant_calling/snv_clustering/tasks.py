'''
Experimental code for clustering allelic counts from multiple samples using a Categorical mixture model fit with VB.

Created on Nov 3, 2015

@author: Andrew Roth
'''
import dask.array as da
import h5py
import numpy as np
import pandas as pd
import yaml

from .batch import BatchVariationalInference
from .svi import StochaticVariationalInference

nucleotides = ['A', 'C', 'G', 'T']

def build_input_file(in_file, out_file, num_dims=4):
    chunk_size = int(1e5)

    hdf_store = pd.HDFStore(in_file, mode='r')
    
    tables = hdf_store.keys()
    
    hdf_store.close()

    with h5py.File(out_file, mode='w') as f:
        for table_name in tables:
            d = f.create_dataset(table_name, compression='gzip', dtype='i8', maxshape=(None, num_dims), shape=(0, num_dims))
            
            for chunk_index, chunk in enumerate(pd.read_hdf(in_file, chunksize=chunk_size, key=table_name, iterator=True)):
                if chunk_index == 0:
                    beg = 0
                    end = chunk.shape[0]
    
                else:
                    beg = end
                    end = beg + chunk.shape[0]
                
                d.resize((end, num_dims))
                
                d[beg:end] = chunk.set_index(['chrom', 'coord']).values

def cluster_snv_counts(
    in_file,
    lower_bound_file=None,
    num_clusters=200,
    params_file=None,
    seed=None):
    
    with h5py.File(in_file, 'r') as f:
        samples = f.keys()
        
        x = []
        
        for key in samples:
            x.append(da.from_array(f[key], chunks=(int(1e4), 4)))
        
        x = da.stack(x, axis=1)
        
        if seed is not None:
            np.random.seed(seed)
        
        gamma_prior = np.ones((x.shape[1], x.shape[2]))
    
        kappa_prior = np.ones(num_clusters)
        
        model_1 = StochaticVariationalInference(gamma_prior, kappa_prior, x, batch_size=1000)
        
        model_1.fit(num_iters=100)
        
        model_2 = BatchVariationalInference(gamma_prior, kappa_prior, x, gamma=model_1.gamma, kappa=model_1.kappa)
        
        model_2.fit()
    
    if lower_bound_file is not None:
        with open(lower_bound_file, 'w') as out_fh:
            yaml.dump({'converged' : model_2.converged, 'lower_bound' : float(model_2.lower_bound[-1])}, out_fh)
    
    if params_file is not None:
        with h5py.File(params_file) as out_fh:
            out_fh.create_dataset('lower_bound', data=model_2.lower_bound)
            
            out_fh.create_dataset('gamma', compression='gzip', data=model_2.gamma)
            
            out_fh.create_dataset('kappa', compression='gzip', data=model_2.kappa)
            
            dset = out_fh.create_dataset('resp', compression='gzip', dtype='f8', shape=(model_2.N, model_2.K))
            
            da.store(model_2.resp, dset)
            

