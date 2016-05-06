from pypeliner.workflow import Workflow


import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.io.download as download
import biowrappers.components.io.download.tasks as download_tasks

import pypeliner

def create_setup_reference_dbs_workflow(config):
    
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
    
    if 'delly' in config:
        workflow.subworkflow(
            name='delly_exclude',
            func=download.create_download_workflow,
            args=(
                config['delly']['exclude_url'],
                pypeliner.managed.OutputFile(config['delly']['exclude_file']),
            )
        )

    if 'destruct' in config:
        import destruct.create_ref_data
        workflow.transform(
            name='destruct_create_ref_data',
            func=destruct.create_ref_data.create_ref_data,
            args=(
                config['destruct']['config'],
                config['destruct']['ref_data_dir'],
            ),
        )

    if 'remixt' in config:
        import remixt.ref_data
        workflow.transform(
            name='remixt_create_ref_data',
            func=remixt.ref_data.create_ref_data,
            args=(
                config['remixt']['config'],
                config['remixt']['ref_data_dir'],
            ),
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

    if 'titan' in config:
        workflow.subworkflow(
            name='gc_wig',
            func=create_gc_wig_file,
            args=(
                config['titan']['config'],
                pypeliner.managed.InputFile(config['ref_genome']['local_path']),
                pypeliner.managed.OutputFile(config['titan']['config']['gc_wig']),
            )
        )

        workflow.subworkflow(
            name='mappability_wig',
            func=create_mappability_wig_file,
            args=(
                config['titan']['config'],
                pypeliner.managed.OutputFile(config['titan']['config']['mappability_wig']),
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
        func=download.create_download_workflow,
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

