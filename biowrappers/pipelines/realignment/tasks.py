'''
Created on Feb 4, 2016

@author: Andrew Roth
'''

import pysam

from biowrappers.components.utils import flatten_input

def write_header_file(in_files, out_file, seq_info):
    
    bam = pysam.AlignmentFile(flatten_input(in_files)[0], mode='r')
    
    header = bam.header.copy()
    
    bam.close()
    
    for x in header['PG'][0]['CL'].split('\t'): 
        if ':' in x:
            key, value = x.split(':')
            
            header['PG'][0][key] = value
            
    header['PG'] = [
        {
            'ID' : 'bwa',
            'VN' : header['PG'][0]['VN'],
            'CL' : 'bwa aln; bwa sampe'
        }
    ]
    
    for entry in header['SQ']:
        for key in seq_info:
            entry[key] = seq_info[key]
    
    pysam.AlignmentFile(out_file, header=header, mode='wh').close()