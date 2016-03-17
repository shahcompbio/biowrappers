import os
import pypeliner
from pypeliner.workflow import Workflow

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.utils as utils

import remixt.workflow

import tasks


def create_remixt_workflow(
    normal_seqdata_file,
    tumour_seqdata_files,
    config,
    out_file,
    raw_data_dir,
    somatic_breakpoint_file=None,
    ref_data_dir=None,
):
    if somatic_breakpoint_file is None:
        raise ValueError('somatic breakpoints required')

    if ref_data_dir is None:
        raise ValueError('ref data directory required')

    results_files = os.path.join(raw_data_dir, 'results', 'sample_{sample_id}.h5')
    selected_files = os.path.join(raw_data_dir, 'selected', 'sample_{sample_id}.h5')
    utils.make_parent_directory(results_files)
    utils.make_parent_directory(selected_files)

    segment_filename = os.path.join(raw_data_dir, 'segment.tsv')

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('sample_id'),
        value=tumour_seqdata_files.keys(),
    )

    workflow.transform(
        name='create_segments',
        ctx={'mem': 4},
        func=remixt.analysis.segment.create_segments,
        args=(
            pypeliner.managed.OutputFile(segment_filename),
            config,
            ref_data_dir,
        ),
        kwargs={
            'breakpoint_filename': pypeliner.managed.InputFile(somatic_breakpoint_file),
        },
    )

    workflow.subworkflow(
        name='remixt',
        func=remixt.workflow.create_remixt_seqdata_workflow,
        args=(
            pypeliner.managed.InputFile(segment_filename),
            pypeliner.managed.InputFile(somatic_breakpoint_file),
            pypeliner.managed.InputFile('tumour_seqdata', 'sample_id', fnames=tumour_seqdata_files),
            pypeliner.managed.InputFile(normal_seqdata_file),
            pypeliner.managed.OutputFile('results', 'sample_id', template=results_files, axes_origin=[]),
            raw_data_dir,
            config,
            ref_data_dir,
        ),
    )

    workflow.transform(
        name='select_solution',
        func=tasks.select_solution,
        axes=('sample_id',),
        args=(
            pypeliner.managed.OutputFile('selected', 'sample_id', template=selected_files),
            pypeliner.managed.InputFile('results', 'sample_id', template=results_files),
        )
    )

    workflow.transform(
        name='merge_results',
        ctx={'mem': 8},
        func=hdf5_tasks.merge_hdf5,
        args=(
            pypeliner.managed.InputFile('selected', 'sample_id', template=selected_files),
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'table_names': '/sample_{}',
        },
    )
    
    return workflow

