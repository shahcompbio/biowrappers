from biowrappers.components.utils import flatten_input

import os
import pypeliner


def merge(in_files, out_file, index_file=None):
    os.environ['MALLOC_ARENA_MAX'] = '4'

    cmd = [
        'picard',
        '-XX:ParallelGCThreads=1',
        '-Xmx8g',
        'MergeSamFiles',
        'OUTPUT={0}'.format(out_file),
        'VALIDATION_STRINGENCY=LENIENT',
    ]

    for file_name in flatten_input(in_files):
        cmd.append('INPUT={0}'.format(file_name))

    pypeliner.commandline.execute(*cmd)

    if index_file is not None:
        cmd = ['samtools', 'index', out_file, index_file]

    pypeliner.commandline.execute(*cmd)


def mark_duplicates(in_files, out_file, metrics_file, index_file=None):
    os.environ['MALLOC_ARENA_MAX'] = '4'

    cmd = [
        'picard',
        '-XX:ParallelGCThreads=1',
        '-Xmx8g',
        'MarkDuplicates',
        'OUTPUT={0}'.format(out_file),
        'METRICS_FILE={0}'.format(metrics_file),
        'VALIDATION_STRINGENCY=LENIENT',
    ]

    for file_name in flatten_input(in_files):
        cmd.append('INPUT={0}'.format(file_name))

    pypeliner.commandline.execute(*cmd)

    if index_file is not None:
        cmd = ['samtools', 'index', out_file, index_file]

    pypeliner.commandline.execute(*cmd)
