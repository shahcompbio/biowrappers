import os
import biowrappers.components.io.download.tasks
import biowrappers.components.io.vcf.tasks

import pypeliner


def download_cosmic(config, out_file, raw_dir):
    for file_type in ('coding', 'non_coding'):
        biowrappers.components.io.download.tasks.download_from_sftp(
            os.path.join(raw_dir, file_type + '.vcf.gz')
            config['host'],
            config['remote_paths'][file_type],
            config['user_name'],
            config['password'])

        pypeliner.commandline.execute(
            'gunzip',
            os.path.join(raw_dir, file_type + '.vcf.gz'))

        biowrappers.components.io.vcf.tasks.finalise_vcf(
            os.path.join(raw_dir, file_type + '.vcf'),
            os.path.join(raw_dir, file_type + '.vcf.gz'))

    biowrappers.components.io.vcf.tasks.concatenate_vcf([
            os.path.join(raw_dir, 'coding.vcf.gz'),
            os.path.join(raw_dir, 'non_coding.vcf.gz')],
        out_file)
