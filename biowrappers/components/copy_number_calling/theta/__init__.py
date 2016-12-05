import os
import pypeliner
from pypeliner.workflow import Workflow

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.utils as utils
import biowrappers.components.io.download

import tasks


def create_theta_workflow(
    seqdata_files,
    config,
    out_file,
    raw_data_dir,
    somatic_breakpoint_file=None,
    normal_id=None,
    **kwargs
):
    if normal_id is None:
        raise ValueError('Theta requires normal sample')

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

    workflow.transform(
        name='run_theta',
        axes=('sample_id',),
        ctx={'mem': 20},
        func=tasks.run_theta,
        args=(
            pypeliner.managed.OutputFile('results', 'sample_id', template=results_files),
            pypeliner.managed.InputFile(normal_seqdata_file),
            pypeliner.managed.InputFile('tumour_seqdata', 'sample_id', fnames=tumour_seqdata_files),
            config,
            pypeliner.managed.TempSpace('work', cleanup=None),
        ),
        kwargs={
            'breakpoints_filename': somatic_breakpoint_file,
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


def create_setup_theta_workflow(config, databases, **kwargs):
    mappability_dir = os.path.realpath(os.path.join(os.path.dirname(config['mappability_template']), os.pardir))
    map_extract_log = os.path.join(mappability_dir, 'mappability_extract.log')
    chromosomes_dir = os.path.dirname(config['chromosome_template'])

    utils.make_directory(mappability_dir)
    utils.make_directory(chromosomes_dir)
    
    workflow = Workflow()

    workflow.subworkflow(
        name='download_mappability', 
        func=biowrappers.components.io.download.create_download_workflow, 
        args=(
            config['mappability_url'],
            pypeliner.managed.TempOutputFile('mappability.tar.gz'),
        )
    )
    
    workflow.commandline(
        name='extract_mappability',
        args=(
            'tar', '-xzvf', pypeliner.managed.TempInputFile('mappability.tar.gz'),
            '-C', mappability_dir,
            '>', pypeliner.managed.OutputFile(map_extract_log),
        ),
    )
    
    for chromosome in config['chromosomes']:
        workflow.subworkflow(
            name='download_chromosome_{}'.format(chromosome),
            func=biowrappers.components.io.download.create_download_workflow, 
            args=(
                config['chromosome_url_template'].format(chromosome),
                pypeliner.managed.TempOutputFile('chromosome_{}.fa.gz'.format(chromosome)),
            )
        )
        
        workflow.commandline(
            name='extract_chromosome_{}'.format(chromosome),
            args=(
                'gunzip', '-c', pypeliner.managed.TempInputFile('chromosome_{}.fa.gz'.format(chromosome)),
                '>', pypeliner.managed.OutputFile(config['chromosome_template'].format(chromosome)),
            ),
        )

    return workflow

