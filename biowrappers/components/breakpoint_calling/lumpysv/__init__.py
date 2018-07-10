import pypeliner
import pypeliner.managed as mgd
from pypeliner.workflow import Workflow

import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.utils as utils
import tasks


def lumpysv_pipeline(
    normal_bam_file,
    tumour_bam_files,
    out_file,
    raw_data_dir,
    normal_sample_id='normal',
):
    bam_files = tumour_bam_files
    bam_files[normal_sample_id] = normal_bam_file

    workflow = pypeliner.workflow.Workflow(default_ctx={'mem': 4})

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=bam_files.keys(),
    )

    workflow.commandline(
        name='extract_discordants',
        axes=('sample_id',),
        args=(
            'samtools', 'view', '-b', '-F', '1294',
            mgd.InputFile('bam', 'sample_id', fnames=bam_files),
            '>',
            mgd.TempOutputFile('discordants.unsorted.bam', 'sample_id'),
        ),
    )

    workflow.commandline(
        name='extract_splitters',
        axes=('sample_id',),
        args=(
            'samtools', 'view', '-h',
            mgd.InputFile('bam', 'sample_id', fnames=bam_files),
            '|',
            'lumpy_extractSplitReads_BwaMem', '-i', 'stdin',
            '|',
            'samtools', 'view', '-Sb', '-',
            '>',
            mgd.TempOutputFile('splitters.unsorted.bam', 'sample_id'),
        ),
    )

    workflow.commandline(
        name='sort_discordants',
        axes=('sample_id',),
        args=(
            'samtools', 'sort',
            mgd.TempInputFile('discordants.unsorted.bam', 'sample_id'),
            '-o',
            mgd.TempOutputFile('discordants.bam', 'sample_id'),
        ),
    )

    workflow.commandline(
        name='sort_splitters',
        axes=('sample_id',),
        args=(
            'samtools', 'sort',
            mgd.TempInputFile('splitters.unsorted.bam', 'sample_id'),
            '-o',
            mgd.TempOutputFile('splitters.bam', 'sample_id'),
        ),
    )

    workflow.transform(
        name='run_lumpyexpress',
        func=tasks.run_lumpyexpress,
        args=(
            mgd.InputFile('bam', 'sample_id', fnames=bam_files),
            mgd.TempInputFile('splitters.bam', 'sample_id'),
            mgd.TempInputFile('discordants.bam', 'sample_id'),
            mgd.TempOutputFile('results.vcf'),
        ),
    )

    workflow.transform(
        name='vcf_to_bcf',
        func=tasks.vcf_to_bcf,
        args=(
            mgd.TempInputFile('results.vcf'),
            mgd.TempOutputFile('results.bcf'),
        ),
    )

    workflow.transform(
        name='convert_bcf',
        func=tasks.convert_bcf,
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        args=(
            mgd.TempInputFile('results.bcf'),
            mgd.OutputFile(out_file),
        ),
        kwargs={
            'control_id': control_id,
        },
    )

    return workflow
