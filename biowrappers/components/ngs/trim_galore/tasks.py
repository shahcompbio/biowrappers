import os
import pypeliner
import shutil


def trim(
        in_fastq_file_1,
        in_fastq_file_2,
        out_fastq_file_1,
        out_fastq_file_2,
        tmp_dir,
        min_length=20,
        out_dir=None,
        run_fastqc=False):

    def get_tmp_file(in_file, read_id, tmp_dir):
        tmp_file = os.path.basename(in_file).split('.')[0]

        tmp_file = tmp_file + '_val_{0}.fq.gz'.format(read_id)

        tmp_file = os.path.join(tmp_dir, tmp_file)

        return tmp_file

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    cmd = [
        'trim_galore',
        '--paired',
        '-o', tmp_dir,
        '--length', min_length,
    ]

    if run_fastqc:
        cmd.append('--fastqc')

    cmd.extend([in_fastq_file_1, in_fastq_file_2])

    pypeliner.commandline.execute(*cmd)

    tmp_fastq_file_1 = get_tmp_file(in_fastq_file_1, 1, tmp_dir)

    shutil.move(tmp_fastq_file_1, out_fastq_file_1)

    tmp_fastq_file_2 = get_tmp_file(in_fastq_file_2, 2, tmp_dir)

    shutil.move(tmp_fastq_file_2, out_fastq_file_2)

    if out_dir is not None:
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        shutil.copytree(tmp_dir, out_dir)


def trim_wrapper(in_files, *args, **kwargs):
    trim(in_files[1], in_files[2], *args, **kwargs)
