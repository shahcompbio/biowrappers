import os
import shutil
import errno
import pypeliner.commandline


def run_optitype(reads_1_fastq, reads_2_fastq, hla_type_file, temp_space):
    try:
        shutil.rmtree(temp_space)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

    pypeliner.commandline.execute(
        'OptiTypePipeline.py', '-i', reads_1_fastq, reads_2_fastq,
        '--dna', '-v', '-o', temp_space)

    results_subdir = os.listdir(temp_space)

    assert len(results_subdir) == 1

    date_time_subdir = results_subdir[0]

    results_file = os.path.join(temp_space, date_time_subdir, date_time_subdir + '_result.tsv')
    os.rename(results_file, hla_type_file)

