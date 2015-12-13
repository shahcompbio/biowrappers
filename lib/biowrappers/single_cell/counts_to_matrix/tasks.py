'''
Created on Oct 29, 2015

@author: andrew
'''
from __future__ import division

import gzip
import pandas as pd

def convert_counts_to_input_matrix(in_file, out_file, col_id_fields, data_type, filters, row_id_fields, min_depth=0, p_value_threshold=1e-6):
    print in_file, out_file, col_id_fields, data_type, filters, row_id_fields
    
    df = pd.read_csv(in_file, compression='gzip', converters={'chrom':str}, sep='\t')
    
    if filters is not None:
        for key, value in filters:
            df = df[df[key] == value]
    
    df['col_id'] = df.apply(lambda row: ':'.join([str(row[field]) for field in col_id_fields]), axis=1)
    
    df['row_id'] = df.apply(lambda row: ':'.join([str(row[field]) for field in row_id_fields]), axis=1)
    
    if 'snv' in data_type:
        df['ref_present'] = df['ref_p_value'] < p_value_threshold
        
        df['alt_present'] = df['alt_p_value'] < p_value_threshold
        
        if data_type == 'binary_snv':
            df['genotype'] = df.apply(_get_binary_snv_genotype, axis=1)
            
        elif data_type == 'snv':
            df['genotype'] = df.apply(_get_snv_genotype, axis=1)
            
        values_col = 'genotype'
        
    elif data_type == 'alt_freq':
        df['depth'] = df['ref_counts'] + df['alt_counts']
    
        df['alt_freq'] = df['alt_counts'] / df['depth']
    
        df.loc[df['depth'] == 0, 'alt_freq'] = pd.np.nan
        
        values_col = 'alt_freq'

    if min_depth > 0:
        df['depth'] = df['ref_counts'] + df['alt_counts']
    
        df.loc[df['depth'] < min_depth, values_col] = pd.np.nan
    
    df = df.pivot(index='row_id', columns='col_id', values=values_col)
    
    df = df.dropna(axis=0, how='all')
    
    df = df.dropna(axis=1, how='all')
    
    with gzip.GzipFile(out_file, 'w') as fh:
        df.to_csv(fh, index_label='cell_id', sep='\t')

def _get_binary_snv_genotype(row):
    if row['alt_present']:
        return 1
    
    elif row['ref_present']:
        return 0
    
    else:
        return pd.np.nan

def _get_snv_genotype(row):
    if row['ref_present'] and row['alt_present']:
        return 1
    
    elif row['ref_present']:
        return 0
    
    elif row['alt_present']:
        return 2
    
    else:
        return pd.np.nan