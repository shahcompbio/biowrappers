'''
Created on Nov 2, 2015

@author: Andrew Roth
'''

import pandas as pd
import pypeliner
import os
import vcf

import biowrappers.components.variant_calling.snpeff.parser

from single_cell.utils import csvutils

def run_snpeff(db, data_dir, in_vcf_file, out_file, classic_mode=True, docker_config={}):

    os.environ['MALLOC_ARENA_MAX'] = '2'

    cmd = [
        'snpEff',
        '-noStats',
        '-noLog',
        '-Xms2g',
        '-Xmx5g',
        '-hgvs1LetterAa',
        '-dataDir',
        data_dir,
    ]

    if classic_mode:
        cmd.append('-classic')

    cmd.extend([
        db,
        in_vcf_file,
        '>',
        out_file
    ])

    pypeliner.commandline.execute(*cmd, **docker_config)


def convert_vcf_to_table(in_file, out_file, table_name, classic_mode=True):
    data = []

    if classic_mode:
        parser = biowrappers.components.variant_calling.snpeff.parser.ClassicSnpEffParser(in_file)

    else:
        parser = biowrappers.components.variant_calling.snpeff.parser.SnpEffParser(in_file)

    for row in parser:
        data.append(row)

    data = pd.DataFrame(data)

    csvutils.write_dataframe_to_csv_and_yaml(data, out_file, data.dtypes.to_dict(),
                                             write_header=True)