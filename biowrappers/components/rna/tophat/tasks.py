from biowrappers.components.utils import make_parent_directory

import pypeliner
import os
import shutil


def align(
    fastq_file_1,
    fastq_file_2,
    ref_genome_index,
    transcriptome_index,
    out_file,
    tmp_dir,
    log_dir=None,
):
    cmd = [
        'tophat',
        '--transcriptome-index', transcriptome_index,
        '--no-coverage-search',
        '-o', tmp_dir,
        ref_genome_index,
        fastq_file_1,
        fastq_file_2
    ]
    pypeliner.commandline.execute(*cmd)
    tmp_out_file = os.path.join(tmp_dir, 'accepted_hits.bam')
    shutil.move(tmp_out_file, out_file)
    if log_dir is not None:
        shutil.copytree(tmp_dir, log_dir)


def build_genome_index(fasta_file, index_prefix):
    make_parent_directory(index_prefix)
    cmd = [
        'bowtie2-build',
        fasta_file,
        index_prefix.replace('.tmp', ''),
    ]
    pypeliner.commandline.execute(*cmd)
    open(index_prefix, 'w').close()


def build_transcriptome_index(
    ref_genome_index_prefix,
    transcript_gtf_file,
    transcriptome_index_prefix,
):
    make_parent_directory(transcriptome_index_prefix)
    cmd = [
        'tophat',
        '-G', transcript_gtf_file,
        '--transcriptome-index', transcriptome_index_prefix.replace('.tmp', ''),
        ref_genome_index_prefix,
    ]
    pypeliner.commandline.execute(*cmd)
    open(transcriptome_index_prefix, 'w').close()
