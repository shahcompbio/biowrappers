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
    bam_files,
    raw_data_dir,
    results_file,
    normal_id=None,
    somatic_breakpoint_file=None,
    patient_config=None,
):
    sample_ids = bam_files.keys()
    
    tumour_ids = bam_files.keys()
    if normal_id is not None:
        tumour_ids.remove(normal_id)

    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('tumour_id'),
        value=tumour_ids,
    )

    seq_data_template = os.path.join(raw_data_dir, 'seqdata', 'sample_{sample_id}.h5')

    if somatic_breakpoint_file is not None:
        somatic_breakpoint_file = pypeliner.managed.InputFile(somatic_breakpoint_file)

    workflow.subworkflow(
        name='extract_seqdata_workflow',
        axes=('sample_id',),
        func=remixt.workflow.create_extract_seqdata_workflow,
        args=(
            pypeliner.managed.InputFile('bam', 'sample_id', fnames=bam_files),
            pypeliner.managed.OutputFile('seqdata', 'sample_id', template=seq_data_template),
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
        assert 'sample_specific' not in remixt_config
        remixt_config.update(patient_config)

        workflow.subworkflow(
            name='remixt',
            func=biowrappers.components.copy_number_calling.remixt.create_remixt_workflow,
            args=(
                pypeliner.managed.InputFile('seqdata', 'sample_id', template=seq_data_template),
                remixt_config,
                pypeliner.managed.OutputFile(remixt_results_filename),
                remixt_raw_data,
            ),
            kwargs={
                'somatic_breakpoint_file': somatic_breakpoint_file,
                'ref_data_dir': config['remixt']['ref_data_dir'],
                'normal_id': normal_id,
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
                pypeliner.managed.InputFile('seqdata', 'sample_id', template=seq_data_template),
                config['titan']['config'],
                pypeliner.managed.OutputFile(titan_results_filename),
                titan_raw_data,
            ),
            kwargs={
                'somatic_breakpoint_file': somatic_breakpoint_file,
                'normal_id': normal_id,
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
                pypeliner.managed.InputFile('seqdata', 'sample_id', template=seq_data_template),
                config['clonehd']['config'],
                pypeliner.managed.OutputFile(clonehd_results_filename),
                clonehd_raw_data,
            ),
            kwargs={
                'somatic_breakpoint_file': somatic_breakpoint_file,
                'normal_id': normal_id,
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

