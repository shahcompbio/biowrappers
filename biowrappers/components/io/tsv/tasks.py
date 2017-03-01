'''
Created on Nov 2, 2015

@author: Andrew Roth
'''

import csv
import gzip
import pandas as pd


def concatenate_tables(in_files, out_file):

    writer = None

    with gzip.GzipFile(out_file, 'w') as out_fh:
        for key in sorted(in_files):
            with gzip.GzipFile(in_files[key], 'r') as in_fh:
                reader = csv.DictReader(in_fh, delimiter='\t')

                if writer is None:
                    writer = csv.DictWriter(out_fh, reader.fieldnames, delimiter='\t')

                    writer.writeheader()

                for row in reader:
                    writer.writerow(row)


def concatenate_tables_hdf5(in_files, out_file, table_name='table'):

    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    for key in sorted(in_files):
        df = pd.read_csv(in_files[key], compression='gzip', sep='\t')

        hdf_store.append(table_name, df)

    hdf_store.close()
