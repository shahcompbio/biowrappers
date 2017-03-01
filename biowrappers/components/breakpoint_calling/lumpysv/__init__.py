import pypeliner
import pypeliner.managed as mgd
from pypeliner.workflow import Workflow

import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.utils as utils
import tasks


def lumpysv_pipeline(
    normal_bam_file,
    tumour_bam_files,
    ref_genome_fasta_file,
    out_file,
    raw_data_dir,
):
    workflow = Workflow()

    workflow.commandline(
        name='extract_unmapped_fastq',
        axes=('sample_id',),
        ctx={'mem': 64, 'num_retry': 2, 'mem_retry_factor': 2},
        args=(
            'samtools', 'view',
            mgd.InputFile('bam', 'sample_id', fnames=bam_filenames),
            '-b', '20',
            '>',
            mgd.TempOutputFile('unmapped_fastq', 'sample_id')
        ),
    )

    workflow.commandline(
        name='bwasw_unmapped',
        axes=('sample_id',),
        ctx={'mem': 64, 'num_retry': 2, 'mem_retry_factor': 2},
        args=(
            'bwa', 'bwasw', '-H', '-t', '20',
            ref_genome_fasta_file,
            mgd.TempInputFile('unmapped_fastq', 'sample_id'),
            '|',
            'samtools', 'view', '-Sb', '-',
            '>',
            mgd.TempOutputFile('split_unsort_bam', 'sample_id')
        ),
    )

                *'bwa bwasw -H -t 20 {genome} {unmapped_fastq} \
            | samtools view -Sb - > {split_unsort_bam}' \


            *'samtools view {paired_bam} \
            | {scripts_directory}/split_unmapped_to_fasta.pl -b 20 \
            > {unmapped_fastq}' \
