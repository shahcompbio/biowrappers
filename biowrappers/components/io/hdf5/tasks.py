'''
Created on Nov 4, 2015

@author: Andrew Roth
'''
import gzip
import pandas as pd

from biowrappers.components.utils import flatten_input

def concatenate_tables(in_files, out_file, non_numeric_as_category=True):
    in_files = flatten_input(in_files)
    
    if non_numeric_as_category:
        col_categories = _get_column_categories(in_files)
    
    else:
        min_itemsize = _get_min_itemsize(in_files)
    
    out_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')
    
    for file_name in in_files:
        in_store = pd.HDFStore(file_name, 'r')
        
        for table_name in in_store.keys():
            df = in_store[table_name]
            
            non_numeric_cols = _get_non_numeric_columns(df)
            
            if non_numeric_as_category:
                for col in non_numeric_cols:
                    df[col] = df[col].astype('category', categories=col_categories[table_name][col])
                
                out_store.append(table_name, df)
            
            else:
                for col in non_numeric_cols:
                    df[col] = df[col].astype(str)
                
                out_store.append(table_name, df, min_itemsize=min_itemsize[table_name])
        
        in_store.close()
    
    out_store.close()

def _get_min_itemsize(file_list):
    
    min_sizes = {}
    
    for file_name in file_list:
        hdf_store = pd.HDFStore(file_name, 'r')
        
        for table_name in hdf_store.keys():
            if table_name not in min_sizes:
                min_sizes[table_name] = {}
            
            df = hdf_store[table_name]
            
            if df.empty:
                continue
            
            for col in _get_non_numeric_columns(df):
                df[col] = df[col].astype(str)
                
                size = max(8, df[col].str.len().max())
                
                if (col not in min_sizes[table_name]) or (size > min_sizes[table_name][col]):
                    min_sizes[table_name][col] = size
        
        hdf_store.close()
        
    return min_sizes

def _get_column_categories(file_list):
    '''
    Find the union set of categories for each column across tables.
    '''
  
    categories = {}
    
    for file_name in file_list:
        hdf_store = pd.HDFStore(file_name, 'r')
        
        for table_name in hdf_store.keys():
            if table_name not in categories:
                categories[table_name] = {}
            
            df = hdf_store[table_name]
            
            if df.empty:
                continue
            
            for col in _get_non_numeric_columns(df):
                df[col] = df[col].astype('category')
                
                if col not in categories[table_name]:
                    categories[table_name][col] = set()
                
                categories[table_name][col].update(set(df[col].cat.categories))
                
        hdf_store.close()
        
    return categories

def _get_non_numeric_columns(df):
    return df.select_dtypes(exclude=[pd.np.number, ]).columns

def convert_hdf5_to_tsv(in_file, key, out_file, compress=False, index=False):
    df = pd.read_hdf(in_file, key)
    
    if compress:
        f_open = gzip.open
    
    else:
        f_open = open
    
    with f_open(out_file, 'w') as fh:
        df.to_csv(fh, index=index, sep='\t')
