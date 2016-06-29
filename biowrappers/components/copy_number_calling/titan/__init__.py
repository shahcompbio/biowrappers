import os
import pypeliner
from pypeliner.workflow import Workflow

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.utils as utils
import biowrappers.components.io.download

import tasks


def create_titan_workflow(
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

    workflow.transform(
        name='prepare_normal_data',
        ctx={'mem': 16, 'num_retry' : 3, 'mem_retry_increment' : 4},
        func=tasks.prepare_normal_data,
        args=(
            pypeliner.managed.InputFile(normal_seqdata_file),
            pypeliner.managed.TempOutputFile('normal.wig'),
            pypeliner.managed.TempOutputFile('het_positions.tsv'),
            config,
        ),
    )

    workflow.transform(
        name='prepare_tumour_data',
        axes=('sample_id',),
        ctx={'mem': 16, 'num_retry' : 3, 'mem_retry_increment' : 4},
        func=tasks.prepare_tumour_data,
        args=(
            pypeliner.managed.InputFile('tumour_seqdata', 'sample_id', fnames=tumour_seqdata_files),
            pypeliner.managed.TempInputFile('het_positions.tsv'),
            pypeliner.managed.TempOutputFile('tumour.wig', 'sample_id'),
            pypeliner.managed.TempOutputFile('tumour_alleles.tsv', 'sample_id'),
            config,
        ),
    )

    workflow.transform(
        name='create_intialization_parameters',
        axes=('sample_id',),
        ctx={'mem': 4, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=tasks.create_intialization_parameters,
        ret=pypeliner.managed.TempOutputObj('init_params', 'sample_id', 'init_param_id'),
        args=(config,),
    )

    workflow.transform(
        name='run_titan',
        axes=('sample_id', 'init_param_id'),
        ctx={'mem': 16, 'num_retry' : 3, 'mem_retry_increment' : 4},
        func=tasks.run_titan,
        args=(
            pypeliner.managed.TempInputObj('init_params', 'sample_id', 'init_param_id'),
            pypeliner.managed.TempInputFile('normal.wig'),
            pypeliner.managed.TempInputFile('tumour.wig', 'sample_id'),
            pypeliner.managed.TempInputFile('tumour_alleles.tsv', 'sample_id'),
            pypeliner.managed.TempOutputFile('cn.tsv', 'sample_id', 'init_param_id'),
            pypeliner.managed.TempOutputFile('params.tsv', 'sample_id', 'init_param_id'),
            config,
        ),
    )

    if somatic_breakpoint_file is not None:
        somatic_breakpoint_file = pypeliner.managed.InputFile(somatic_breakpoint_file)

    workflow.transform(
        name='select_solution',
        axes=('sample_id',),
        ctx={'mem': 4, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=tasks.select_solution,
        args=(
            pypeliner.managed.TempInputObj('init_params', 'sample_id', 'init_param_id'),
            pypeliner.managed.TempInputFile('cn.tsv', 'sample_id', 'init_param_id'),
            pypeliner.managed.TempInputFile('params.tsv', 'sample_id', 'init_param_id'),
            pypeliner.managed.OutputFile('results', 'sample_id', template=results_files),
            pypeliner.managed.OutputFile(os.path.join(raw_data_dir, 'output', '{sample_id}_cn_loci.tsv'), 'sample_id'),
            pypeliner.managed.OutputFile(os.path.join(raw_data_dir, 'output', '{sample_id}_cn_segments.tsv'), 'sample_id'),
            pypeliner.managed.OutputFile(os.path.join(raw_data_dir, 'output', '{sample_id}_cn_igv.tsv'), 'sample_id'),
            pypeliner.managed.OutputFile(os.path.join(raw_data_dir, 'output', '{sample_id}_params.tsv'), 'sample_id'),
            config,
            pypeliner.managed.Template('{sample_id}', 'sample_id'),
        ),
        kwargs={
            'breakpoints_filename': somatic_breakpoint_file,
        },
    )

    workflow.transform(
        name='merge_results',
        ctx={'mem': 8, 'num_retry': 3, 'mem_retry_increment': 2},
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


def create_gc_wig_file(config, genome_file, out_file):
    workflow = Workflow()

    workflow.commandline(
        name='create_gc',
        ctx={'mem': 4},
        args=(
            'gcCounter',
            '-w', config['window_size'],
            pypeliner.managed.InputFile(genome_file),
            '>',
            pypeliner.managed.OutputFile(out_file),
        ),
    )

    return workflow


def create_mappability_wig_file(config, out_file):
    workflow = Workflow()
    
    workflow.subworkflow(
        name='download_mappability_bigwig',
        func=biowrappers.components.io.download.create_download_workflow,
        args=(
            config['mappability_url'],
            pypeliner.managed.TempOutputFile('mappability_bigwig'),
        )
    )

    workflow.commandline(
        name='convert_mappability_to_wig',
        ctx={'mem': 4},
        args=(
            'mapCounter',
            '-w', config['window_size'],
            pypeliner.managed.TempInputFile('mappability_bigwig'),
            '>',
            pypeliner.managed.OutputFile(out_file),
        ),
    )

    return workflow


def create_setup_titan_workflow(config, databases, **kwargs):
    workflow = Workflow()

    workflow.subworkflow(
        name='gc_wig',
        func=create_gc_wig_file,
        args=(
            config,
            pypeliner.managed.InputFile(databases['ref_genome']['local_path']),
            pypeliner.managed.OutputFile(config['gc_wig']),
        )
    )

    workflow.subworkflow(
        name='mappability_wig',
        func=create_mappability_wig_file,
        args=(
            config,
            pypeliner.managed.OutputFile(config['mappability_wig']),
        )
    )

    return workflow

