import os
import pypeliner
from pypeliner.workflow import Workflow

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.utils as utils

import tasks


def create_battenberg_single_workflow(
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
        ctx={'mem': 20},
        func=tasks.prepare_data,
        args=(
            pypeliner.managed.InputFile(normal_seqdata_file),
            pypeliner.managed.InputFile(tumour_seqdata_file),
            pypeliner.managed.TempOutputFile('allele_counts.tar.gz'),
            config,
        ),
    )

    if somatic_breakpoint_file is not None:
        somatic_breakpoint_file = pypeliner.managed.InputFile(somatic_breakpoint_file)

    workflow.transform(
        name='run_battenberg',
        ctx={'mem': 8},
        func=tasks.run_battenberg,
        args=(
            pypeliner.managed.TempInputFile('allele_counts.tar.gz'),
            pypeliner.managed.OutputFile(results_file),
            pypeliner.managed.TempSpace('run_battenberg_temp', cleanup=None),
            config
        ),
        kwargs={
            'somatic_breakpoint_file': somatic_breakpoint_file,
        },
    )

    return workflow


def create_battenberg_workflow(
    seqdata_files,
    config,
    out_file,
    raw_data_dir,
    somatic_breakpoint_file=None,
    normal_id=None,
    **kwargs
):
    if normal_id is None:
        raise ValueError('cloneHD requires normal sample')

    normal_seqdata_file = seqdata_files[normal_id]
    tumour_seqdata_files = seqdata_files.copy()
    del tumour_seqdata_files[normal_id]

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
        name='run_battenberg',
        axes=('sample_id',),
        func=create_battenberg_single_workflow,
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


def create_setup_battenberg_workflow(config, databases, **kwargs):
    return Workflow()
