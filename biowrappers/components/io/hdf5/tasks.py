'''
Created on Nov 4, 2015

@author: Andrew Roth
'''
from collections import defaultdict

import gzip
import pandas as pd
import re

from biowrappers.components.utils import flatten_input
from pandas.api.types import CategoricalDtype

def concatenate_tables(
    in_files,
    out_file,
    drop_duplicates=False,
    in_memory=True,
    non_numeric_as_category=True
):
    in_files = flatten_input(in_files)

    # Only support drop duplicatess in memory
    if drop_duplicates or in_memory:
        _concatenate_tables_in_memory(
            in_files,
            out_file,
            drop_duplicates=drop_duplicates,
            non_numeric_as_category=non_numeric_as_category
        )

    else:
        _concatenate_tables_on_disk(
            in_files,
            out_file,
            non_numeric_as_category=non_numeric_as_category
        )


def _concatenate_tables_in_memory(
    in_files,
    out_file,
    drop_duplicates=False,
    non_numeric_as_category=True
):

    tables = defaultdict(list)

    for file_name in in_files:
        in_store = pd.HDFStore(file_name, 'r')

        for table_name in _iter_table_names(in_store):
            df = in_store[table_name]
            tables[table_name].append(df)

        in_store.close()

    for table_name in tables:
        # TODO: check columns / types are equivalent
        non_empty_tables = list(filter(lambda a: not a.empty, tables[table_name]))

        if len(non_empty_tables) == 0:
            tables[table_name] = tables[table_name][0]
            continue

        else:
            tables[table_name] = pd.concat(non_empty_tables)

        if drop_duplicates:
            tables[table_name] = tables[table_name].drop_duplicates()

        if non_numeric_as_category:
            non_numeric_cols = _get_non_numeric_columns(tables[table_name])

            for col in non_numeric_cols:
                if pd.api.types.infer_dtype(tables[table_name][col]) == 'unicode':
                    tables[table_name][col] = tables[table_name][col].astype(str)
                tables[table_name][col] = tables[table_name][col].astype('category')

    out_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    for table_name in tables:
        if tables[table_name].empty:
            out_store.put(table_name, tables[table_name])
        else:
            out_store.put(table_name, tables[table_name], format='table')

    out_store.close()


def _concatenate_tables_on_disk(in_files, out_file, non_numeric_as_category=True):
    if non_numeric_as_category:
        col_categories = _get_column_categories(in_files)

    else:
        min_itemsize = _get_min_itemsize(in_files)

    out_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    table_columns = defaultdict(set)

    for file_name in in_files:
        in_store = pd.HDFStore(file_name, 'r')

        for table_name in _iter_table_names(in_store):
            df = in_store[table_name]

            table_columns[table_name].update(df.columns)

            if df.empty:
                continue

            non_numeric_cols = _get_non_numeric_columns(df)

            if non_numeric_as_category:
                for col in non_numeric_cols:
                    if df[col].dtype.name == 'category':
                        df[col] = df[col].cat.set_categories(col_categories[table_name][col])

                    else:
                        df[col] = df[col].astype(CategoricalDtype(categories=col_categories[table_name][col]))

                out_store.append(table_name, df, format='table')

            else:
                for col in non_numeric_cols:
                    df[col] = df[col].astype(str)

                out_store.append(table_name, df, min_itemsize=min_itemsize[table_name], format='table')

        for table_name, columns in table_columns.items():
            out_store.append(table_name, pd.DataFrame(columns=columns), format='table')

        in_store.close()

    out_store.close()


def _get_min_itemsize(file_list):
    '''
    Get the minimum string size for all columns in a list of tables from a list of HDFStores.
    '''

    min_sizes = {}

    for file_name in file_list:
        hdf_store = pd.HDFStore(file_name, 'r')

        for table_name in _iter_table_names(hdf_store):
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

        for table_name in _iter_table_names(hdf_store):
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
    '''
    Find the set of non-numeric (int, float, complex) columns in a table.
    '''

    return df.select_dtypes(exclude=[pd.np.number, ]).columns


def _iter_table_names(store):
    '''
    Returns an iterator over all non-metadata tables in Pandas HDFStore.
    '''

    meta_data_tables = _get_meta_data_tables(store)

    for table_name in store.keys():
        if table_name in meta_data_tables:
            continue

        yield table_name


def _get_meta_data_tables(store):
    '''
    Find all tables in and HDFStore which are pandas meta-data tables.
    '''

    meta_tables = set()

    for table_name in store.keys():
        for meta_table_name in store.keys():
            if (re.search('.*/meta/.*/meta$', meta_table_name) is not None) and meta_table_name.startswith(table_name):
                meta_tables.add(meta_table_name)

    return meta_tables


def convert_hdf5_to_tsv(in_file, key, out_file, compress=False, index=False):
    '''
    Convert and pandas HDF5 table to tsv format.
    '''

    df = pd.read_hdf(in_file, key)

    if compress:
        f_open = gzip.open

    else:
        f_open = open

    with f_open(out_file, 'w') as fh:
        df.to_csv(fh, index=index, sep='\t')


def merge_hdf5(in_files, out_file, table_names='{}'):
    '''
    Merge pandas HDF5 tables
    '''

    out_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    for file_key, file_name in in_files.items():
        in_store = pd.HDFStore(file_name, 'r')

        # Compatability with dictionary keyed by single or multiple values
        if not isinstance(file_key, tuple):
            file_key = (file_key,)

        for table_name in _iter_table_names(in_store):
            df = in_store[table_name]

            # Workaround: currently cannot store empty dataframe in table format
            format = 'table'
            if len(df.index) == 0:
                format = None

            out_store.put(table_names.format(*file_key) + '/' + table_name, df, format=format)
