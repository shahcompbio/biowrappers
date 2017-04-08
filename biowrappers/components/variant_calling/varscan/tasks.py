'''
Created on 8 Apr 2017

@author: Andrew Roth
'''
import os
import pypeliner.commandline as cli
import shutil


def run_somatic(normal_file, tumour_file, out_file, tmp_dir):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_prefix = os.path.join(tmp_dir, 'varscan')

    cmd = [
        'varscan',
        'somatic',
        normal_file,
        tumour_file,
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

    cli.execute('bcftools', 'concat', snp_file, indel_file, '-a', '-o', out_file, '-O', 'v')

    shutil.rmtree(tmp_dir)
