import gzip
import itertools
import os
import pipes
import pypeliner
import pysam
import shutil
import subprocess

from biowrappers.components.utils import flatten_input


def concatenate(in_files, out_file):
    """ Concatenate FASTQ files.
    """
    with gzip.GzipFile(out_file, 'w') as out_fh:
        for in_file in flatten_input(in_files):
            with gzip.GzipFile(in_file, 'r') as in_fh:
                out_fh.write(in_fh.read())


def is_sanger(file_name, num_reads=10000):
    quals = set()
    with pysam.FastxFile(file_name) as reader:
        for i, record in enumerate(reader):
            if i >= num_reads:
                break
            quals.update(set([ord(x) for x in record.quality]))
    return min(quals) < 64


def convert_qualities_to_sanger(in_file, out_file):
    if is_sanger(in_file):
        os.link(in_file, out_file)
        #shutil.copyfile(in_file, out_file)
    else:
        cmd = [
            'gzip', '-cd', in_file,
            '|',
            'sed', '-e',
            "'4~4y/@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\\\]^_`abcdefghi/!\"#$%&'\\''()*+,-.\\/0123456789:;<=>?@ABCDEFGHIJ/'",
            '|',
            'gzip', '-cf',
            '>',
            out_file
        ]
        #pypeliner.commandline.execute(*cmd)
        subprocess.check_call(' '.join(cmd), shell=True)

def split_fastq(in_filename, out_filenames, num_reads_per_file):
    """ Split a fastq file.
    """

    if in_filename.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    with opener(in_filename, 'r') as in_file:
        file_number = 0
        out_file = None
        out_file_read_count = None
        try:
            for name, seq, comment, qual in itertools.izip_longest(*[in_file] * 4):
                if out_file is None or out_file_read_count == num_reads_per_file:
                    if out_file is not None:
                        out_file.close()
                    out_file = open(out_filenames[file_number], 'w')
                    out_file_read_count = 0
                    file_number += 1
                out_file.write(name)
                out_file.write(seq)
                out_file.write(comment)
                out_file.write(qual)
                out_file_read_count += 1
        finally:
            if out_file is not None:
                out_file.close()
