'''
Created on Feb 3, 2016

@author: Andrew Roth
'''
import os
import pypeliner


def run_aln(in_fastq_file, ref_genome_fasta_file, out_sai_file):
    pypeliner.commandline.execute(
        'bwa',
        'aln',
        ref_genome_fasta_file,
        in_fastq_file,
        '>',
        out_sai_file
    )
    assert os.path.getsize(out_sai_file) > 0


def run_sampe(
        in_fastq_file_1,
        in_fastq_file_2,
        in_sai_file_1,
        in_sai_file_2,
        ref_genome_fasta_file,
        out_file,
        out_file_format='bam',
        num_compression_threads=0,
        read_group_info=None):

    cmd = [
        'bwa',
        'sampe'
    ]

    if read_group_info is not None:
        read_group_str = ['@RG', 'ID:{0}'.format(read_group_info['ID'])]

        for key, value in sorted(read_group_info.items()):
            if key == 'ID':
                continue

            read_group_str.append(':'.join((key, value)))

        read_group_str = '\t'.join(read_group_str)

        cmd.extend(['-r', read_group_str])

    cmd.extend([
        ref_genome_fasta_file,
        in_sai_file_1,
        in_sai_file_2,
        in_fastq_file_1,
        in_fastq_file_2
    ])

    cmd.extend(['|', 'samtools', 'view', '-@', num_compression_threads])

    if out_file_format == 'bam':
        cmd.extend(
            [
                '-b'
            ]
        )

    elif out_file_format == 'cram':
        cmd.extend(
            [
                '-C',
                '-T', ref_genome_fasta_file
            ]
        )

    else:
        raise ValueError('{0} is not a valid file format.'.format(out_file_format))

    cmd.extend(['>', out_file])

    pypeliner.commandline.execute(*cmd)
