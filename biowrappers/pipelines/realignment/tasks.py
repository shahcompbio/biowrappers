'''
Created on Feb 4, 2016

@author: Andrew Roth
'''

import pysam

from biowrappers.components.utils import flatten_input


def get_read_group_config(file_name):
    bam = pysam.AlignmentFile(file_name, mode='rb', check_sq=False)

    assert len(bam.header['RG']) == 1

    config = bam.header['RG'][0].copy()

    bam.close()

    return config


def write_header_file(in_files, out_file, seq_info):

    bam = pysam.AlignmentFile(flatten_input(in_files)[0], mode='r', check_sq=False)

    header = bam.header.copy()

    bam.close()

    for x in header['PG'][0]['CL'].split('\t'):
        if ':' in x:
            key, value = x.split(':')

            header['PG'][0][key] = value

    header['PG'] = [
        {
            'ID': 'bwa',
            'VN': header['PG'][0]['VN'],
            'CL': 'bwa aln; bwa sampe'
        }
    ]

    for entry in header['SQ']:
        for key in seq_info:
            entry[key] = seq_info[key]

    pysam.AlignmentFile(out_file, header=header, mode='wh').close()
