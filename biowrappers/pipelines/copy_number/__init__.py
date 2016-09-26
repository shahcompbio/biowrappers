from pypeliner.workflow import Workflow

import os
import pypeliner

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
from biowrappers.components.utils import make_parent_directory

import biowrappers.components.copy_number_calling.remixt
import biowrappers.components.copy_number_calling.titan
import biowrappers.components.copy_number_calling.clonehd

import remixt.workflow


def call_and_annotate_pipeline(
    config,
    normal_bam_file,
    tumour_bam_files,
    raw_data_dir,
    results_file,
    normal_id='normal',
    somatic_breakpoint_file=None,
    ploidy_config=None,
):
    sample_ids = tumour_bam_files.keys()

    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('tumour_id'),
        value=sample_ids,
    )

    tumour_seq_data_template = os.path.join(raw_data_dir, 'seqdata', 'sample_{tumour_id}.h5')
    normal_seq_data_filename = os.path.join(raw_data_dir, 'seqdata', 'sample_{}.h5'.format(normal_id))

    if somatic_breakpoint_file is not None:
        somatic_breakpoint_file = pypeliner.managed.InputFile(somatic_breakpoint_file)

    workflow.subworkflow(
        name='extract_seqdata_workflow_normal',
        func=remixt.workflow.create_extract_seqdata_workflow,
        args=(
            pypeliner.managed.InputFile(normal_bam_file),
            pypeliner.managed.OutputFile(normal_seq_data_filename),
            config['remixt'].get('extract_seqdata', {}),
            config['remixt']['ref_data_dir'],
        ),
    )

    workflow.subworkflow(
        name='extract_seqdata_workflow_tumour',
        axes=('tumour_id',),
        func=remixt.workflow.create_extract_seqdata_workflow,
        args=(
            pypeliner.managed.InputFile('bam', 'tumour_id', fnames=tumour_bam_files),
            pypeliner.managed.OutputFile('seqdata', 'tumour_id', template=tumour_seq_data_template),
            config['remixt'].get('extract_seqdata', {}),
            config['remixt']['ref_data_dir'],
        ),
    )

    merge_inputs = {}
    
    if 'remixt' in config:
        remixt_raw_data = os.path.join(raw_data_dir, 'remixt')
        remixt_results_filename = os.path.join(remixt_raw_data, 'results.h5')
        make_parent_directory(remixt_results_filename)

        remixt_config = config['remixt']['config']
        if ploidy_config is not None:
            if 'sample_specific' not in remixt_config:
                remixt_config['sample_specific'] = {}
            remixt_config['sample_specific'].update(ploidy_config)

        workflow.subworkflow(
            name='remixt',
            func=biowrappers.components.copy_number_calling.remixt.create_remixt_workflow,
            args=(
                pypeliner.managed.InputFile(normal_seq_data_filename),
                pypeliner.managed.InputFile('seqdata', 'tumour_id', template=tumour_seq_data_template),
                remixt_config,
                pypeliner.managed.OutputFile(remixt_results_filename),
                remixt_raw_data,
            ),
            kwargs={
                'somatic_breakpoint_file': somatic_breakpoint_file,
                'ref_data_dir': config['remixt']['ref_data_dir'],
            },
        )

        merge_inputs['/copy_number/remixt'] = pypeliner.managed.InputFile(remixt_results_filename)

    if 'titan' in config:
        titan_raw_data = os.path.join(raw_data_dir, 'titan')
        titan_results_filename = os.path.join(titan_raw_data, 'results.h5')
        make_parent_directory(titan_results_filename)

        workflow.subworkflow(
            name='titan',
            func=biowrappers.components.copy_number_calling.titan.create_titan_workflow,
            args=(
                pypeliner.managed.InputFile(normal_seq_data_filename),
                pypeliner.managed.InputFile('seqdata', 'tumour_id', template=tumour_seq_data_template),
                config['titan']['config'],
                pypeliner.managed.OutputFile(titan_results_filename),
                titan_raw_data,
            ),
            kwargs={
                'somatic_breakpoint_file': somatic_breakpoint_file,
            },
        )

        merge_inputs['/copy_number/titan'] = pypeliner.managed.InputFile(titan_results_filename)

    if 'clonehd' in config:
        clonehd_raw_data = os.path.join(raw_data_dir, 'clonehd')
        clonehd_results_filename = os.path.join(clonehd_raw_data, 'results.h5')
        make_parent_directory(clonehd_results_filename)

        workflow.subworkflow(
            name='clonehd',
            func=biowrappers.components.copy_number_calling.clonehd.create_clonehd_workflow,
            args=(
                pypeliner.managed.InputFile(normal_seq_data_filename),
                pypeliner.managed.InputFile('seqdata', 'tumour_id', template=tumour_seq_data_template),
                config['clonehd']['config'],
                pypeliner.managed.OutputFile(clonehd_results_filename),
                clonehd_raw_data,
            ),
            kwargs={
                'somatic_breakpoint_file': somatic_breakpoint_file,
            },
        )

        merge_inputs['/copy_number/clonehd'] = pypeliner.managed.InputFile(clonehd_results_filename)

    workflow.transform(
        name='merge_results',
        ctx={'mem': 8},
        func=hdf5_tasks.merge_hdf5,
        args=(
            merge_inputs,
            pypeliner.managed.OutputFile(results_file),
        ),
    )

    return workflow

