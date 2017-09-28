'''
Created on Nov 2, 2015

@author: Andrew Roth
'''

import pandas as pd
import pypeliner
import os
import vcf

import biowrappers.components.variant_calling.snpeff.parser


def run_snpeff(db, in_vcf_file, out_file, classic_mode=True):

    os.environ['MALLOC_ARENA_MAX'] = '2'

    cmd = [
        'snpEff',
        '-noStats',
        '-noLog',
        '-Xms2g',
        '-Xmx5g',
        '-hgvs1LetterAa',
    ]

    if classic_mode:
        cmd.append('-classic')

    cmd.extend([
        db,
        in_vcf_file,
        '>',
        out_file
    ])

    pypeliner.commandline.execute(*cmd)


def convert_vcf_to_table(in_file, out_file, table_name, classic_mode=True):
    data = []

    if classic_mode:
        parser = biowrappers.components.variant_calling.snpeff.parser.ClassicSnpEffParser(in_file)

    else:
        parser = biowrappers.components.variant_calling.snpeff.parser.SnpEffParser(in_file)

    for row in parser:
        data.append(row)

    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    hdf_store[table_name] = pd.DataFrame(data)

    hdf_store.close()
