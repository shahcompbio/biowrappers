'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
import os
import pypeliner
import time


def run_mutect(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        cosmic_vcf_file,
        dbsnp_vcf_file,
        region,
        out_file):

    cmd = [
        'gatk',
        '-T', 'MuTect2',
        '-I:normal', normal_bam_file,
        '-I:tumor', tumour_bam_file,
        '-R', ref_genome_fasta_file,
        '--cosmic', cosmic_vcf_file,
        '--dbsnp', dbsnp_vcf_file,
        '-o', out_file,
        '-L', region,
        '>', '/dev/null'  # Avoid dumping mutect tsv file to stdout
    ]

    pypeliner.commandline.execute(*cmd)

    idx_file_name = out_file + '.idx'

    # Guard against slow file system
    time.sleep(10)
    
    if os.path.getsize(out_file) == 0:
        os.unlink(out_file)

        raise Exception('{} is empty. Removing.'.format(out_file))

    os.unlink(idx_file_name)
