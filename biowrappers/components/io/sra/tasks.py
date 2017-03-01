import os
import pypeliner.commandline


def sra_fastq_dump_paired(sra_filename, fastq_1_filename, fastq_2_filename, temp_directory):
    """ Dump paired fastq files from SRA format to 2 fastq files.

    Args:
        sra_filename (str): path to sra filename
        fastq_1_filename (str): path to output fastq filename for end 1
        fastq_2_filename (str): path to output fastq filename for end 2
        temp_directory (str): temporary working directory
    """

    pypeliner.commandline.execute(
        'fastq-dump',
        '--split-files',
        '--origfmt',
        '--outdir', temp_directory,
        sra_filename,
    )

    assert sra_filename.endswith('.sra')

    temp_fq_1_filenames = os.path.join(temp_directory, os.path.basename(sra_filename)[:-4] + '_1.fastq')
    temp_fq_2_filenames = os.path.join(temp_directory, os.path.basename(sra_filename)[:-4] + '_2.fastq')

    os.rename(temp_fq_1_filenames, fastq_1_filename)
    os.rename(temp_fq_2_filenames, fastq_2_filename)
