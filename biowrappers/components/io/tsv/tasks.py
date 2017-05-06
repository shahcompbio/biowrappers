'''
Created on Nov 2, 2015

@author: Andrew Roth
'''
from biowrappers.components.utils import flatten_input
from pandas.errors import EmptyDataError

import csv
import gzip
import pandas as pd


def concatenate_tables(in_files, out_file, ignore_empty_files=False, use_gzip=True):

    if use_gzip:
        open_func = gzip.GzipFile
    
    else:
        open_func = open

    write_header = True
    
    with open_func(out_file, 'w') as out_fh:
        for file_name in flatten_input(in_files):
            try:
                df = pd.read_csv(file_name, sep='\t')
            
            except EmptyDataError as e:
                if ignore_empty_files:
                    continue

                else:
                    raise e

            if df.empty:
                continue
 
            df.to_csv(out_fh, header=write_header, index=False, sep='\t')
            
            write_header = False


def concatenate_tables_hdf5(in_files, out_file, table_name='table'):

    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    for key in sorted(in_files):
        df = pd.read_csv(in_files[key], compression='gzip', sep='\t')

        hdf_store.append(table_name, df)

    hdf_store.close()
