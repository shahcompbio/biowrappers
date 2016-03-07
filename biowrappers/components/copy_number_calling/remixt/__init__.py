import os
import pypeliner
import pypeliner.managed as mgd
from pypeliner.workflow import Workflow

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.utils as utils

import remixt.workflow


def remixt_pipeline(
    normal_bam_file,
    tumour_bam_files,
    somatic_breakpoint_filename,
    config,
    ref_data_dir,
    out_file,
    raw_data_dir,
    normal_sample_id='normal',
):
    results_files = os.path.join(raw_data_dir, 'results', 'sample_{sample_id}.h5')
    utils.make_parent_directory(results_files)

    segment_filename = os.path.join(raw_data_dir, 'segment.tsv')

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('sample_id'),
        value=tumour_bam_files.keys(),
    )

    workflow.transform(
        name='create_segments',
        ctx={'mem': 4},
        func=remixt.analysis.segment.create_segments,
        args=(
            mgd.OutputFile(segment_filename),
            config,
        ),
        kwargs={
            'breakpoint_filename': mgd.InputFile(somatic_breakpoint_filename),
        },
    )

    workflow.subworkflow(
        name='remixt_pipeline',
        func=remixt.workflow.create_remixt_pipeline,
        args=(
            mgd.InputFile(segment_filename),
            mgd.InputFile(somatic_breakpoint_filename),
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumour_bam_files),
            mgd.InputFile(normal_bam_file),
            normal_sample_id,
            mgd.OutputFile(results_files, 'sample_id', axes_origin=[]),
            raw_data_dir,
            mgd.TempInputObj('remixt_config'),
        ),
    )

    merge_inputs = {}
    for sample_id in tumour_bam_files.keys():
        merge_inputs['/sample_{sample_id}'.format(sample_id)] = pypeliner.managed.InputFile(results_files.format(sample_id=sample_id))

    workflow.transform(
        name='merge_results',
        ctx={'mem': 8},
        func=hdf5_tasks.merge_hdf5,
        args=(
            merge_inputs,
            pypeliner.managed.OutputFile(out_file),
        )
    )
    
    return workflow

