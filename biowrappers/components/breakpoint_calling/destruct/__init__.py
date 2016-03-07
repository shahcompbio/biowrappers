import os
import pypeliner
import pypeliner.managed as mgd
from pypeliner.workflow import Workflow

import biowrappers.components.utils as utils
import tasks

import destruct.workflow


def destruct_pipeline(
    normal_bam_file,
    tumour_bam_files,
    config,
    ref_data_dir,
    out_file,
    raw_data_dir,
    normal_sample_id='normal',
):
    bam_files = tumour_bam_files
    bam_files[normal_sample_id] = normal_bam_file

    utils.make_directory(os.path.join(raw_data_dir, 'raw'))
    breakpoint_file = os.path.join(raw_data_dir, 'raw', 'breakpoint.tsv')
    breakpoint_library_file = os.path.join(raw_data_dir, 'raw', 'breakpoint_library.tsv')

    utils.make_directory(os.path.join(raw_data_dir, 'somatic'))
    somatic_breakpoint_file = os.path.join(raw_data_dir, 'somatic', 'breakpoint.tsv')
    somatic_breakpoint_library_file = os.path.join(raw_data_dir, 'somatic', 'breakpoint_library.tsv')

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('sample_id'),
        value=bam_files.keys(),
    )

    workflow.subworkflow(
        name='run_destruct',
        func=destruct.workflow.create_destruct_workflow,
        args=(
            pypeliner.managed.InputFile('bam', 'sample_id', fnames=bam_files),
            pypeliner.managed.OutputFile(breakpoint_file),
            pypeliner.managed.OutputFile(breakpoint_library_file),
            config,
            ref_data_dir,
        ),
    )

    workflow.transform(
        name='filter_breakpoints',
        ctx={'mem': 4},
        func=tasks.filter_breakpoints,
        args=(
            pypeliner.managed.InputFile(breakpoint_file),
            pypeliner.managed.InputFile(breakpoint_library_file),
            normal_sample_id,
            pypeliner.managed.OutputFile(somatic_breakpoint_file),
            pypeliner.managed.OutputFile(somatic_breakpoint_library_file),
        ),
    )

    workflow.transform(
        name='write_store',
        func=tasks.write_store,
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        args=(
            pypeliner.managed.InputFile(somatic_breakpoint_file),
            pypeliner.managed.InputFile(somatic_breakpoint_library_file),
            mgd.OutputFile(out_file),
        ),
    )
    
    return workflow

