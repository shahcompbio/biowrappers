from pypeliner.workflow import Workflow

import os
import pypeliner

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
from biowrappers.components.utils import make_parent_directory

import remixt.workflow


def call_and_annotate_pipeline(
    config,
    normal_bam_file,
    tumour_bam_files,
    somatic_breakpoint_file,
    raw_data_dir,
    results_file,
):
    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('tumour_sample_id'),
        value=tumour_bam_files.keys(),
    )

    merge_inputs = {}
    
    if 'remixt' in config:
        remixt_raw_data = os.path.join(raw_data_dir, 'remixt')
        remixt_results_filename = os.path.join(remixt_raw_data, 'results.h5')
        make_parent_directory(remixt_results_filename)

        workflow.subworkflow(
            name='remixt',
            func=remixt.remixt_pipeline,
            args=(
                pypeliner.managed.InputFile(normal_bam_file),
                pypeliner.managed.InputFile('tumour_bams', 'tumour_sample_id', fnames=tumour_bam_files),
                pypeliner.managed.InputFile(somatic_breakpoint_file),
                config['remixt']['config'],
                config['remixt']['ref_data_dir'],
                pypeliner.managed.OutputFile(remixt_results_filename),
                remixt_raw_data,
            ),
        )

        merge_inputs['/copy_number/remixt'] = pypeliner.managed.InputFile(remixt_results_filename)

    workflow.transform(
        name='merge_results',
        ctx={'mem': 8},
        func=hdf5_tasks.merge_hdf5,
        args=(
            merge_inputs,
            pypeliner.managed.OutputFile(results_file),
        )
    )

    return workflow

