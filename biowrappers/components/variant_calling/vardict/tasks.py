'''
Created on Nov 1, 2015

@author: Andrew Roth
'''

import pipes
import pypeliner
import subprocess
import vcf


def run_single_sample_vardict(
        bam_file,
        ref_genome_fasta_file,
        region,
        out_file,
        java=False,
        min_allele_frequency=0.01,
        remove_duplicate_reads=False,
        sample_name=None):
    try:
        region = region.split('-')
        if ':' in region[0]:
            region = region[0].split(':') + region[1:]
    except:
        pass

    if java:
        prog = 'vardict-java'
    else:
        prog = 'vardict'
    cmd = [
        prog,
        '-b', bam_file,
        '-f', min_allele_frequency,
        '-G', ref_genome_fasta_file,
        '-R', '{0}:{1}-{2}'.format(*region),
    ]
    if java:
        cmd.extend(['-th', 1])
    if remove_duplicate_reads:
        cmd.append('-t')
    cmd.extend(['|', 'teststrandbias.R', '|', 'var2vcf_valid.pl'])
    if sample_name is not None:
        cmd.extend(['-N', sample_name])
    cmd.extend(['>', out_file])
    pypeliner.commandline.execute(*cmd)


def run_paired_sample_vardict(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        region,
        out_file,
        java=False,
        min_allele_frequency=0.01,
        remove_duplicate_reads=False,
        sample_names=None):
    """ Run VarDict in paired mode

    :param sample_names: dictionary with sample names. Should have keys `normal` and `tumour`.
    """
    try:
        region = region.split('-')
        if ':' in region[0]:
            region = region[0].split(':') + region[1:]
    except:
        pass

    if java:
        prog = 'vardict-java'
    else:
        prog = 'vardict'
    cmd = [
        prog,
        '-b', pipes.quote('{0}|{1}'.format(tumour_bam_file, normal_bam_file)),
        '-f', min_allele_frequency,
        '-G', ref_genome_fasta_file,
        '-R', '{0}:{1}-{2}'.format(*region),
    ]
    if java:
        cmd.extend(['-th', 1])
    if remove_duplicate_reads:
        cmd.append('-t')
    cmd.extend(['|', 'testsomatic.R', '|', 'var2vcf_paired.pl', '-f', min_allele_frequency])
    if sample_names is not None:
        cmd.extend(['-N', pipes.quote('{tumour}|{normal}'.format(**sample_names))])
    cmd.extend(['>', out_file])
    cmd_str = ' '.join([str(x) for x in cmd])
    subprocess.check_call(cmd_str, shell=True)


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
