'''
Created on 8 Apr 2017

@author: Andrew Roth
'''
import os
import pypeliner.commandline as cli
import shutil

from biowrappers.components.ngs.samtools.tasks import mpileup


def run_pileup2snp(in_file, out_file):
    cmd = [
        'varscan',
        'pileup2snp',
        in_file,
        '|',
        'gzip',
        '>',
        out_file
    ]

    cli.execute(*cmd)


def run_somatic(normal_bam, tumour_bam, ref_genome_fasta_file, out_file, region, tmp_dir):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    normal_pileup = os.path.join(tmp_dir, 'normal.pileup')

    tumour_pileup = os.path.join(tmp_dir, 'tumour.pileup')

    mpileup(normal_bam, normal_pileup, ref_genome_fasta_file, region)

    mpileup(tumour_bam, tumour_pileup, ref_genome_fasta_file, region)

    tmp_prefix = os.path.join(tmp_dir, 'varscan')

    cmd = [
        'varscan',
        'somatic',
        normal_pileup,
        tumour_pileup,
        tmp_prefix,
        '--output-vcf'
    ]

    cli.execute(*cmd)

    snp_file = tmp_prefix + '.snp.vcf'

    indel_file = tmp_prefix + '.indel.vcf'

    cli.execute('bgzip', snp_file)

    cli.execute('bgzip', indel_file)

    snp_file = snp_file + '.gz'

    indel_file = indel_file + '.gz'

    cli.execute('tabix', snp_file)

    cli.execute('tabix', indel_file)

    cli.execute('bcftools', 'concat', snp_file, indel_file, '-a', '-o', out_file, '-O', 'z')

    shutil.rmtree(tmp_dir)
