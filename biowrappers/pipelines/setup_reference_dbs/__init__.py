from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.io.download as download
import biowrappers.components.io.download.tasks as download_tasks

def create_init_reference_dbs_workflow(config):
    
    workflow = Workflow()
    
    if 'cosmic' in config:
        workflow.subworkflow(
            name='cosmic',
            func=create_cosmic_download_workflow,
            args=(
                config['cosmic'],
                pypeliner.managed.OutputFile(config['cosmic']['local_path']),
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
    
    if 'ref_genome' in config:
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

    workflow.subworkflow(
        name='delly_exclude', 
        func=download.create_download_workflow, 
        args=(
            config['delly']['exclude_url'],
            pypeliner.managed.OutputFile(config['delly']['exclude_file']),
        )
    )

    return workflow

def create_cosmic_download_workflow(config, out_file):
    
    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('remote_path', 'files'),
        value={
            'coding' : config['remote_paths']['coding'],
            'non_coding' : config['remote_paths']['non_coding']
        },
    )
    
    workflow.transform(
        name='download_files',
        axes=('files',),
        ctx={'local' : True},
        func=download_tasks.download_from_sftp,
        args=(
            pypeliner.managed.TempOutputFile('download.vcf.gz', 'files'),
            config['host'],
            pypeliner.managed.TempInputObj('remote_path', 'files'),
            config['user_name'],
            config['password']
        )
    )
    
    workflow.transform(
        name='merge',
        axes=(),
        ctx={'mem' : 4},
        func=vcf_tasks.concatenate_vcf_fast,
        args=(
            pypeliner.managed.TempInputFile('download.vcf.gz', 'files'),
            pypeliner.managed.TempOutputFile('merged.vcf'),
        ),
    )
    
    workflow.transform(
        name='sort',
        axes=(),
        ctx={'mem' : 4},
        func=vcf_tasks.sort_vcf,
        args=(
            pypeliner.managed.TempInputFile('merged.vcf'),
            pypeliner.managed.TempOutputFile('sorted.vcf')
        )
    )
    
    workflow.subworkflow(
        name='finalise',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('sorted.vcf'),
            pypeliner.managed.OutputFile(out_file)
        ),
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
            pypeliner.managed.InputFile(out_file)
        )
    )
    
    workflow.commandline(
        name='build_fai',
        ctx={'mem' : 6, 'num_retry' : 3, 'mem_retry_increment' : 2},
        args=(
            'samtools',
            'faidx',
            pypeliner.managed.InputFile(out_file)
        )
    )
    
    workflow.commandline(
        name='build_bwa_index', 
        ctx={'mem' : 6, 'num_retry' : 3, 'mem_retry_increment' : 2},
        args=(
            'bwa',
            'index',
            pypeliner.managed.InputFile(out_file)
        )
    )
    
    return workflow
