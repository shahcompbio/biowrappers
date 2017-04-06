import os
import pypeliner
import shutil


def build_index(index_file, transcriptome_fasta_file, kmer_length=31):
    cmd = [
        'kallisto',
        'index',
        '-i', index_file,
        '-k', kmer_length,
        transcriptome_fasta_file,
    ]

    pypeliner.commandline.execute(*cmd)


def quantify(
        fastq_file_1,
        fastq_file_2,
        index_file,
        out_file,
        tmp_dir,
        num_bootstraps=100):

    cmd = [
        'kallisto',
        'quant',
        '-i', index_file,
        '-o', tmp_dir,
        '-b', num_bootstraps,
        fastq_file_1,
        fastq_file_2,
    ]

    pypeliner.commandline.execute(*cmd)

    tmp_out_file = os.path.join(tmp_dir, 'abundance.tsv')

    shutil.move(tmp_out_file, out_file)


def quantify_wrapper(fastq_files, *args, **kwargs):
    quantify(fastq_files[1], fastq_files[2], *args, **kwargs)
