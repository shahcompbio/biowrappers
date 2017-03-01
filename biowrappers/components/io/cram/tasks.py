'''
Created on Feb 4, 2016

@author: Andrew Roth
'''

import pypeliner

from biowrappers.components.utils import flatten_input


def merge(
        in_files,
        reference_genome_fasta_file,
        out_file,
        attach_read_group_from_file_name=False,
        header_file=None,
        compression_level=9,
        num_compression_threads=0):

    cmd = [
        'samtools',
        'merge',
        '-c',
        '-p',
        '-f',
        '-l', compression_level,
        '-@', num_compression_threads,
        '--output-fmt', 'cram',
        '--reference', reference_genome_fasta_file,
    ]

    if attach_read_group_from_file_name:
        cmd.append('-r')

    if header_file is not None:
        cmd.extend(['-h', header_file])

    cmd.append(out_file)

    for file_name in flatten_input(in_files):
        cmd.append(file_name)

    pypeliner.commandline.execute(*cmd)


def sort(
        in_file,
        reference_genome_fasta_file,
        out_file,
        max_mem='2G',
        name_sort=False,
        compression_level=9,
        num_threads=1):

    cmd = [
        'samtools',
        'sort',
        '-l', compression_level,
        '-m', max_mem,
        '-o', out_file,
        '-@', num_threads,
        '--output-fmt', 'cram',
        '--reference', reference_genome_fasta_file
    ]

    if name_sort:
        cmd.append('-n')

    cmd.append(in_file)

    pypeliner.commandline.execute(*cmd)
