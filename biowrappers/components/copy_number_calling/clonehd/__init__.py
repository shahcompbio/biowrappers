import os
import pypeliner
from pypeliner.workflow import Workflow

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.utils as utils

import tasks


def create_clonehd_single_workflow(
    normal_seqdata_file,
    tumour_seqdata_file,
    config,
    results_file,
    somatic_breakpoint_file=None,
    **kwargs
):
    workflow = Workflow()

    workflow.transform(
        name='prepare_data',
        ctx={'mem': 16},
        func=tasks.prepare_data,
        args=(
            pypeliner.managed.InputFile(normal_seqdata_file),
            pypeliner.managed.InputFile(tumour_seqdata_file),
            pypeliner.managed.TempOutputFile('normal.cna.txt'),
            pypeliner.managed.TempOutputFile('tumour.cna.txt'),
            pypeliner.managed.TempOutputFile('tumour.baf.txt'),
            config,
        ),
    )

    workflow.transform(
        name='run_clonehd',
        ctx={'mem': 8},
        func=tasks.run_clonehd,
        args=(
            pypeliner.managed.TempInputFile('normal.cna.txt'),
            pypeliner.managed.TempInputFile('tumour.cna.txt'),
            pypeliner.managed.TempInputFile('tumour.baf.txt'),
            pypeliner.managed.TempOutputFile('tumour.summary.txt'),
            pypeliner.managed.TempOutputFile('cna_subclone', 'subclone'),
            pypeliner.managed.TempOutputFile('bam_subclone', 'subclone'),
            pypeliner.managed.TempSpace('run_clonehd_temp'),
        ),
    )

    if somatic_breakpoint_file is not None:
        somatic_breakpoint_file = pypeliner.managed.InputFile(somatic_breakpoint_file)

    workflow.transform(
        name='report',
        ctx={'mem': 4},
        func=tasks.report,
        args=(
            pypeliner.managed.TempInputFile('tumour.summary.txt'),
            pypeliner.managed.TempInputFile('cna_subclone', 'subclone'),
            pypeliner.managed.TempInputFile('bam_subclone', 'subclone'),
            pypeliner.managed.OutputFile(results_file),
        ),
        kwargs={
            'somatic_breakpoint_file': somatic_breakpoint_file,
        },
    )

    return workflow


def create_clonehd_workflow(
    normal_seqdata_file,
    tumour_seqdata_files,
    config,
    out_file,
    raw_data_dir,
    somatic_breakpoint_file=None,
    **kwargs
):
    results_files = os.path.join(raw_data_dir, 'results', 'sample_{sample_id}.h5')
    utils.make_parent_directory(results_files)

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('sample_id'),
        value=tumour_seqdata_files.keys(),
    )

    if somatic_breakpoint_file is not None:
        somatic_breakpoint_file = pypeliner.managed.InputFile(somatic_breakpoint_file)

    workflow.subworkflow(
        name='run_clonehd',
        axes=('sample_id',),
        func=create_clonehd_single_workflow,
        args=(
            pypeliner.managed.InputFile(normal_seqdata_file),
            pypeliner.managed.InputFile('tumour_seqdata', 'sample_id', fnames=tumour_seqdata_files),
            config,
            pypeliner.managed.OutputFile('results', 'sample_id', template=results_files),
        ),
        kwargs={
            'somatic_breakpoint_file': somatic_breakpoint_file,
        },
    )

    workflow.transform(
        name='merge_results',
        ctx={'mem': 8},
        func=hdf5_tasks.merge_hdf5,
        args=(
            pypeliner.managed.InputFile('results', 'sample_id', template=results_files),
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'table_names': '/sample_{}',
        },
    )

    return workflow


def create_setup_clonehd_workflow(config, databases, **kwargs):
    return Workflow()

