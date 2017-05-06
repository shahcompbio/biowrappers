import os
import pypeliner.commandline
import shutil


def sra_fastq_dump_paired(sra_filename, fastq_1_filename, fastq_2_filename, tmp_dir, clean_tmp_dir=False):
    """ Dump paired fastq files from SRA format to 2 fastq files.

    Args:
        sra_filename (str): path to sra filename
        fastq_1_filename (str): path to output fastq filename for end 1
        fastq_2_filename (str): path to output fastq filename for end 2
        tmp_dir (str): temporary working directory
    """
    if os.path.exists(tmp_dir) and clean_tmp_dir:
        shutil.rmtree(tmp_dir)

    pypeliner.commandline.execute(
        'fastq-dump',
        '--gzip',
        '--split-files',
        '--origfmt',
        '--outdir', tmp_dir,
        sra_filename,
    )

    assert sra_filename.endswith('.sra')

    temp_fq_1_filenames = os.path.join(tmp_dir, os.path.basename(sra_filename)[:-4] + '_1.fastq.gz')
    temp_fq_2_filenames = os.path.join(tmp_dir, os.path.basename(sra_filename)[:-4] + '_2.fastq.gz')

    os.rename(temp_fq_1_filenames, fastq_1_filename)
    os.rename(temp_fq_2_filenames, fastq_2_filename)

    if clean_tmp_dir:
        shutil.rmtree(tmp_dir)
