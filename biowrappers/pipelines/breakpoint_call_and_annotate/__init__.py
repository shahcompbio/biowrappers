from pypeliner.workflow import Workflow

import os
import pypeliner  

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
from biowrappers.components.utils import make_parent_directory

import biowrappers.components.breakpoint_calling.destruct as destruct


def call_and_annotate_pipeline(
    config,
    normal_bam_path,
    tumour_bam_paths,
    raw_data_dir,
    results_file,
):
    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('tumour_sample_id'),
        value=tumour_bam_paths.keys(),
    )

    merge_inputs = {}
    
    if 'destruct' in config:
        destruct_raw_data = os.path.join(raw_data_dir, 'destruct')
        destruct_results_filename = os.path.join(destruct_raw_data, 'results.h5')
        make_parent_directory(destruct_results_filename)

        workflow.subworkflow(
            name='destruct',
            func=destruct.destruct_pipeline,
            args=(
                pypeliner.managed.InputFile(normal_bam_path),
                pypeliner.managed.InputFile('tumour_bams', 'tumour_sample_id', fnames=tumour_bam_paths),
                config['destruct']['config'],
                config['destruct']['ref_data_dir'],
                pypeliner.managed.OutputFile(destruct_results_filename),
                destruct_raw_data,
            ),
        )

        merge_inputs['/breakpoints/destruct'] = pypeliner.managed.InputFile(destruct_results_filename)

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

