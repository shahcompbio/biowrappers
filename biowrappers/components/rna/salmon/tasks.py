from biowrappers.components.utils import make_parent_directory

import os
import pypeliner
import shutil


def build_index(
        index_sentinel_file,
        transcriptome_fasta_file,
        kmer_length=31,
        gencode=False,
        num_threads=1):

    make_parent_directory(index_sentinel_file)

    cmd = [
        'salmon',
        'index',
        '-i', os.path.dirname(index_sentinel_file),
        '-k', kmer_length,
        '-p', num_threads,
        '-t', transcriptome_fasta_file,
    ]

    if gencode is not None:
        cmd.append('--gencode')

    pypeliner.commandline.execute(*cmd)

    open(index_sentinel_file, 'w').close()


def quantify(
        fastq_file_1,
        fastq_file_2,
        index,
        out_file,
        tmp_dir,
        library_type='A',
        num_threads=1):

    cmd = [
        'salmon',
        'quant',
        '-i', index,
        '-1', fastq_file_1,
        '-2', fastq_file_2,
        '-l', library_type,
        '-p', num_threads,
        '-o', tmp_dir,
    ]

    pypeliner.commandline.execute(*cmd)

    tmp_out_file = os.path.join(tmp_dir, 'quant.sf')

    shutil.move(tmp_out_file, out_file)


def quantify_wrapper(fastq_files, *args, **kwargs):
    quantify(fastq_files[1], fastq_files[2], *args, **kwargs)
