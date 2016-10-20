'''
Created on Nov 1, 2015

@author: Andrew Roth
'''

import pypeliner
import vcf


def run_single_sample_vardict(
        bam_file,
        ref_genome_fasta_file,
        region,
        out_file,
        min_allele_frequency=0.01,
        sample_name=None):
    cmd = [
        'vardict',
        '-b', bam_file,
        '-f', min_allele_frequency,
        '-G', ref_genome_fasta_file,
        '-R', '{0}:{1}-{2}'.format(*region),
        '-th', 1,
    ]
    cmd.extend(['|', 'teststrandbias.R', '|', 'var2vcf_valid.pl'])
    if sample_name is not None:
        cmd.extend(['-N', sample_name])
    cmd.extend(['>', out_file])    
    pypeliner.commandline.execute(*cmd)


def run_vardict(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        out_file,
        min_allele_frequency=0.01,
        region=None,
        tumour_sample_name=None):

    cmd = [
        'vardict-java',
        '-b', '{0}|{1}'.format(tumour_bam_file, normal_bam_file),
        '-f', min_allele_frequency,
        '-G', ref_genome_fasta_file,
        '-th', 1
    ]

    if region is not None:
        cmd.extend(['-R', region])

    if tumour_sample_name is not None:
        cmd.extend(['-N', tumour_sample_name])

    cmd.extend(['>', out_file])

    pypeliner.commandline.execute(*cmd)


def run_vardict_test_somatic(in_file, out_file):

    cmd = [
        'testsomatic.R',
        '<', in_file,
        '>', out_file
    ]

    pypeliner.commandline.execute(*cmd)


def run_vardict_var_to_vcf(
        in_file,
        out_file,
        min_allele_frequency=0.01):

    cmd = [
        'var2vcf_paired.pl',
        '-f', min_allele_frequency,
        '<', in_file,
        '>', out_file
    ]

    pypeliner.commandline.execute(*cmd)


def filter_vcf(in_file, out_file, variant_type):

    reader = vcf.Reader(filename=in_file)

    with open(out_file, 'w') as out_fh:
        writer = vcf.Writer(out_fh, reader)

        for record in reader:
            if record.INFO['STATUS'] != 'StrongSomatic':
                continue

            if len(record.FILTER) > 0:
                continue

            if (variant_type == 'indel') and (record.is_indel):
                writer.write_record(record)

            elif (variant_type == 'snv') and (record.is_snp):
                writer.write_record(record)
