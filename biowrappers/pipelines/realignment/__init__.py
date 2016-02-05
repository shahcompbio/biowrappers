from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.bam.tasks as bam_tasks
import biowrappers.components.alignment.bwa.tasks as bwa_tasks

import tasks

def realignment_pipeline(
        config,
        in_file, 
        out_file):

    ref_genome = pypeliner.managed.InputFile(config['ref_genome']['file'])
    
    read_1 = pypeliner.managed.TempFile('read_1', 'split')
    
    read_2 = pypeliner.managed.TempFile('read_2', 'split')
    
    read_1_sai = pypeliner.managed.TempFile('read_1.sai', 'split')
    
    read_2_sai = pypeliner.managed.TempFile('read_2.sai', 'split')
    
    workflow = Workflow()
    
    workflow.transform(
        name='bam_to_fasta',
        axes=(),
        ctx={'mem' : 4, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=bam_tasks.convert_to_fastqs, 
        args=(
            pypeliner.managed.InputFile(in_file),
            {
                1 : read_1.as_output(),
                2 : read_2.as_output(),
            },
            pypeliner.managed.TempSpace('bam_to_fastq'),
        ),
        kwargs={
            'split_size' : config['split_size']
        },
    )
    
    workflow.transform(
        name='aln_read_1',
        axes=('split',),
        ctx={'mem' : 4, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=bwa_tasks.run_aln,
        args=(
            read_1.as_input(),
            ref_genome,
            read_1_sai.as_output(),
        ),
    )
    
    workflow.transform(
        name='aln_read_2',
        axes=('split',),
        ctx={'mem' : 4, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=bwa_tasks.run_aln,
        args=(
            read_2.as_input(),
            ref_genome,
            read_2_sai.as_output(),
        ),
    )
    
    workflow.transform(
        name='sampe', 
        axes=('split',),
        ctx={'mem' : 4, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=bwa_tasks.run_sampe, 
        args=(
            read_1.as_input(),
            read_2.as_input(),
            read_1_sai.as_input(),
            read_2_sai.as_input(),
            ref_genome,
            pypeliner.managed.TempOutputFile('aligned.bam', 'split'),
        ),
        kwargs={
            'read_group_info' : config['read_group']
        },
    )
    
    workflow.transform(
        name='sort',
        axes=('split',),
        ctx={'mem' : 4, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=bam_tasks.sort,
        args=(
            pypeliner.managed.TempInputFile('aligned.bam', 'split'),
            pypeliner.managed.TempOutputFile('sorted.bam', 'split'),
        ),
    )
    
    workflow.transform(
        name='write_header_file',
        axes=(),
        ctx={'local' : True},
        func=tasks.write_header_file,
        args=(
            pypeliner.managed.TempInputFile('sorted.bam', 'split'),
            pypeliner.managed.TempOutputFile('header.sam'),
            config['ref_genome']['header']
        ),
    )
    
    workflow.transform(
        name='merge',
        axes=(),
        ctx={'mem' : 4, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=bam_tasks.merge, 
        args=(
            pypeliner.managed.TempInputFile('sorted.bam', 'split'),
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'header_file' : pypeliner.managed.TempInputFile('header.sam'),
        },
    )
    
    return workflow
