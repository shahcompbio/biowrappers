from biowrappers.components.utils import flatten_input

import pypeliner


def merge(in_files, out_file):
    cmd = [
        'picard',
        'MergeSamFiles',
        'SO=coordinate',
        'USE_THREADING=false',
        'VALIDATION_STRINGENCY=LENIENT',
        'OUTPUT={}'.format(out_file),
    ]
    cmd.extend(['INPUT={}'.format(x) for x in flatten_input(in_files)])
    pypeliner.commandline.execute(*cmd)
