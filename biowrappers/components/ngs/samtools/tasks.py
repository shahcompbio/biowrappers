import gzip
import pandas as pd
import pypeliner
import shutil


def depth(in_file, out_file):
    """ Run `samtools depth`. Output will be gzip compressed.
    """
    cmd = [
        'samtools', 'depth', in_file,
        '|', 'gzip',
        '>', out_file,
    ]
    pypeliner.commandline.execute(*cmd)


def depth_to_bed(in_file, out_file, min_depth=1):
    """ Convert the output of `samtools depth` to a bed format file. Output will be gzip compressed BED file.
    """
    df = pd.read_csv(in_file, compression='gzip', header=None, names=['chrom', 'coord', 'depth'], sep='\t')
    df['beg'] = df['coord'] - 1
    df['end'] = df['coord']
    df = df[df['depth'] >= min_depth]
    df = df[['chrom', 'beg', 'end', 'depth']]
    with gzip.GzipFile(out_file, 'w') as out_fh:
        df.to_csv(out_fh, header=False, index=False, sep='\t')


def faidx(in_file, out_file):
    cmd = [
        'samtools',
        'faidx',
        in_file,
    ]
    pypeliner.commandline.execute(*cmd)
    shutil.move(in_file + '.fai', out_file)
