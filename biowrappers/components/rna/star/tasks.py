''''
Created on 20 Jul 2016

@author: Andrew Roth
'''
from biowrappers.components.utils import make_parent_directory

import os
import pypeliner
import shutil


def align(
        fastq_file_1,
        fastq_file_2,
        ref_genome_dir,
        out_file,
        tmp_dir,
        log_dir=None,
        num_threads=1,
        read_group_info=None,
        unaligned_read_fastq_1=None,
        unaligned_read_fastq_2=None):

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_prefix = os.path.join(tmp_dir, 'alignment')

    cmd = [
        'STAR',
        '--runThreadN', num_threads,
        '--genomeDir', ref_genome_dir,
        '--readFilesIn', fastq_file_1, fastq_file_2,
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--outSAMattributes', 'NH', 'HI', 'NM', 'MD', 'AS', 'nM',
        '--readFilesCommand', 'zcat',
        '--outFileNamePrefix', tmp_prefix,
    ]

    if read_group_info is not None:
        read_group_str = ['ID:{0}'.format(read_group_info['ID']), ]
        for key, value in sorted(read_group_info.items()):
            if key == 'ID':
                continue

            read_group_str.append(':'.join((key, value)))

        read_group_str = '\t'.join(read_group_str)

        cmd.extend(['--outSAMattrRGline', read_group_str])

    if unaligned_read_fastq_1 is not None:
        cmd.extend(['--outReadsUnmapped', 'Fastx'])

    pypeliner.commandline.execute(*cmd)

    tmp_out_file = tmp_prefix + 'Aligned.sortedByCoord.out.bam'

    shutil.move(tmp_out_file, out_file)

    if unaligned_read_fastq_1 is not None:
        tmp_out_file = tmp_prefix + 'Unmapped.out.mate1'

        shutil.move(tmp_out_file, unaligned_read_fastq_1)

    if unaligned_read_fastq_2 is not None:
        tmp_out_file = tmp_prefix + 'Unmapped.out.mate2'

        shutil.move(tmp_out_file, unaligned_read_fastq_2)

    if log_dir is not None:
        if os.path.exists(log_dir):
            shutil.rmtree(log_dir)

        shutil.copytree(tmp_dir, log_dir)


def build_index(
        index_sentinel_file,
        ref_genome_fasta_file,
        transcript_gtf_file,
        overhang=100,
        num_threads=1):

    make_parent_directory(index_sentinel_file)

    cmd = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--runThreadN', num_threads,
        '--genomeDir', os.path.dirname(index_sentinel_file),
        '--genomeFastaFiles', ref_genome_fasta_file,
        '--sjdbGTFfile', transcript_gtf_file,
        '--sjdbOverhang', overhang,
    ]

    pypeliner.commandline.execute(*cmd)

    open(index_sentinel_file, 'w').close()
