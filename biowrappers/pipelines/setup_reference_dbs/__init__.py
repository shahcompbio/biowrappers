from pypeliner.workflow import Workflow


import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.io.download as download
import biowrappers.components.io.download.tasks as download_tasks
import biowrappers.components.copy_number_calling.remixt
import biowrappers.components.copy_number_calling.theta
import biowrappers.components.copy_number_calling.titan
import biowrappers.components.copy_number_calling.clonehd

import tasks

import pypeliner


def create_setup_reference_dbs_workflow(config):
    
    workflow = Workflow()
    
    if 'cosmic' in config:
        workflow.transform(
            name='cosmic',
            func=tasks.download_cosmic,
            args=(
                config['cosmic'],
                pypeliner.managed.OutputFile(config['cosmic']['local_path']),
                pypeliner.managed.TempSpace('cosmic_work', cleanup=None)
            )
        )
    
    if 'dbsnp' in config:
        workflow.subworkflow(
            name='dbsnp',
            func=create_dbsnp_download_workflow,
            args=(
                config['dbsnp'],
                pypeliner.managed.OutputFile(config['dbsnp']['local_path']),
            )
        )
    
    if 'mappability' in config:
        workflow.subworkflow(
            name='mappability', 
            func=download.create_download_workflow, 
            args=(
                config['mappability']['url'],
                pypeliner.managed.OutputFile(config['mappability']['local_path']),
            )
        )
    
    if 'ref_genome' in config and 'url' in config['ref_genome']:
        workflow.subworkflow(
            name='ref_genome', 
            func=create_ref_genome_download_and_index_workflow, 
            args=(
                config['ref_genome'],
                pypeliner.managed.OutputFile(config['ref_genome']['local_path']),
            )
        )
        
    if 'snpeff' in config:
        workflow.commandline(
            name='snpeff', 
            args=(
                'snpEff',
                'download',
                config['snpeff']['db']
            )
        )

    if 'chrom_info' in config:
        workflow.subworkflow(
            name='chrom_info',
            func=download.create_download_workflow,
            args=(
                config['chrom_info']['url'],
                pypeliner.managed.OutputFile(config['chrom_info']['local_path']),
            )
        )

    return workflow


def create_setup_tools_workflow(databases, config):
    
    workflow = Workflow()

    if 'destruct' in config:
        import destruct.create_ref_data
        workflow.transform(
            name='destruct_create_ref_data',
            ctx={'mem': 16},
            func=destruct.create_ref_data.create_ref_data,
            args=(
                config['destruct']['config'],
                config['destruct']['ref_data_dir'],
            ),
        )

    if 'delly' in config:
        workflow.subworkflow(
            name='delly_exclude',
            func=download.create_download_workflow,
            args=(
                config['delly']['exclude_url'],
                pypeliner.managed.OutputFile(config['delly']['exclude_file']),
            )
        )

    if 'remixt' in config:
        workflow.subworkflow(
            name='create_setup_remixt_workflow',
            func=biowrappers.components.copy_number_calling.remixt.create_setup_remixt_workflow,
            args=(
                config['remixt']['config'],
                databases,
            ),
            kwargs={
                'ref_data_dir': config['remixt']['ref_data_dir'],
            },
        )

    if 'titan' in config:
        workflow.subworkflow(
            name='create_setup_titan_workflow',
            func=biowrappers.components.copy_number_calling.titan.create_setup_titan_workflow,
            args=(
                config['titan']['config'],
                databases,
            )
        )

    if 'theta' in config:
        workflow.subworkflow(
            name='create_setup_theta_workflow',
            func=biowrappers.components.copy_number_calling.theta.create_setup_theta_workflow,
            args=(
                config['theta']['config'],
                databases,
            )
        )

    if 'clonehd' in config:
        workflow.subworkflow(
            name='create_setup_clonehd_workflow',
            func=biowrappers.components.copy_number_calling.clonehd.create_setup_clonehd_workflow,
            args=(
                config['clonehd']['config'],
                databases,
            )
        )

    return workflow


def create_dbsnp_download_workflow(config, out_file):
    
    workflow = Workflow()
    
    workflow.subworkflow(
        name='download',
        func=download.create_download_workflow,
        args=(
            config['url'], 
            pypeliner.managed.OutputFile(out_file)
        )
    )
        
    workflow.transform(
        name='index',
        ctx={'mem' : 4},
        func=vcf_tasks.index_vcf,
        args=(
            pypeliner.managed.InputFile(out_file),
            pypeliner.managed.OutputFile(out_file + '.tbi')
        )
    )
    
    return workflow

def create_ref_genome_download_and_index_workflow(config, out_file):

    workflow = Workflow()
    
    workflow.subworkflow(
        name='download',
        func=download.create_download_workflow,
        args=(
            config['url'], 
            pypeliner.managed.OutputFile(out_file)
        )
    )
    
    workflow.commandline(
        name='build_dict',
        ctx={'mem' : 6, 'num_retry' : 3, 'mem_retry_increment' : 2},
        args=(
            'samtools',
            'dict',
            pypeliner.managed.InputFile(out_file),
            '>',
            pypeliner.managed.OutputFile(out_file + '.build_dict.log'),
        )
    )
    
    workflow.commandline(
        name='build_fai',
        ctx={'mem' : 6, 'num_retry' : 3, 'mem_retry_increment' : 2},
        args=(
            'samtools',
            'faidx',
            pypeliner.managed.InputFile(out_file),
            '>',
            pypeliner.managed.OutputFile(out_file + '.build_fai.log'),
        )
    )
    
    workflow.commandline(
        name='build_bwa_index', 
        ctx={'mem' : 6, 'num_retry' : 3, 'mem_retry_increment' : 2},
        args=(
            'bwa',
            'index',
            pypeliner.managed.InputFile(out_file),
            '>',
            pypeliner.managed.OutputFile(out_file + '.build_bwa_index.log'),
        )
    )
    
    return workflow

