from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.bam.tasks as bam_tasks
import biowrappers.components.alignment.bwa.tasks as bwa_tasks

import tasks


def realignment_pipeline(
        config,
        in_file,
        out_file,
        read_group_info=None):

    if read_group_info is None:
        read_group_info = config.get('read_group', {})

    if 'ID' not in read_group_info:
        read_group_info['ID'] = hash(in_file) % int(1e6)

    ref_genome = pypeliner.managed.InputFile(config['ref_genome']['file'])

    read_1 = pypeliner.managed.TempFile('read_1', 'split')

    read_2 = pypeliner.managed.TempFile('read_2', 'split')

    read_1_sai = pypeliner.managed.TempFile('read_1.sai', 'split')

    read_2_sai = pypeliner.managed.TempFile('read_2.sai', 'split')

    read_group_config = pypeliner.managed.TempObj('read_group_config')

    workflow = Workflow()

    if 'read_group' in config:
        workflow.setobj(
            obj=read_group_config.as_output(),
            value=read_group_info,
        )

    else:
        workflow.transform(
            name='get_read_group_config',
            ctx={'local': True},
            func=tasks.get_read_group_config,
            ret=read_group_config.as_output(),
            args=(
                pypeliner.managed.InputFile(in_file),
            )
        )

    workflow.transform(
        name='bam_to_fasta',
        axes=(),
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=bam_tasks.convert_to_fastqs,
        args=(
            pypeliner.managed.InputFile(in_file),
            {
                1: read_1.as_output(),
                2: read_2.as_output(),
            },
            pypeliner.managed.TempSpace('bam_to_fastq'),
        ),
        kwargs={
            'split_size': config['split_size']
        },
    )

    workflow.transform(
        name='aln_read_1',
        axes=('split',),
        ctx={'mem': 6, 'num_retry': 3, 'mem_retry_increment': 2},
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
        ctx={'mem': 6, 'num_retry': 3, 'mem_retry_increment': 2},
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
        ctx={'mem': 6, 'num_retry': 3, 'mem_retry_increment': 2},
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
            'read_group_info': read_group_config.as_input()
        },
    )

    workflow.transform(
        name='sort',
        axes=('split',),
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=bam_tasks.sort,
        args=(
            pypeliner.managed.TempInputFile('aligned.bam', 'split'),
            pypeliner.managed.TempOutputFile('sorted.bam', 'split'),
        ),
    )

    workflow.transform(
        name='write_header_file',
        axes=(),
        ctx={'local': True},
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
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=bam_tasks.merge,
        args=(
            pypeliner.managed.TempInputFile('sorted.bam', 'split'),
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'header_file': pypeliner.managed.TempInputFile('header.sam'),
        },
    )

    return workflow


def realignment_readgroups_pipeline(
        config,
        in_file,
        out_file):

    workflow = Workflow()

    workflow.transform(
        name='get_read_group_configs',
        func=tasks.get_read_group_configs,
        ret=pypeliner.managed.TempOutputObj('read_group_id', 'read_group_config'),
        args=(
            pypeliner.managed.InputFile(in_file),
        )
    )

    workflow.commandline(
        name='create_read_group_bam',
        axes=('read_group_id',),
        args=(
            'samtools', 'view', '-b',
            '-r', pypeliner.managed.InputInstance('read_group_id'),
            pypeliner.managed.InputFile(in_file),
            '>',
            pypeliner.managed.TempOutputFile('read_group_bam', 'read_group_id'),
        )
    )

    workflow.subworkflow(
        name='realignment_pipeline',
        axes=('read_group_id',),
        func=realignment_pipeline,
        args=(
            config,
            pypeliner.managed.TempInputFile('read_group_bam', 'read_group_id'),
            pypeliner.managed.TempOutputFile('realigned_read_group_bam', 'read_group_id'),
        )
        kwargs={
            'read_group_info': pypeliner.managed.TempOutputObj('read_group_config', 'read_group_id'),
        }
    )

    workflow.transform(
        name='merge_and_markdups',
        axes=('read_group_id',),
        ctx={'mem' : 48, 'num_retry' : 3, 'mem_retry_increment' : 16},
        func=bam_tasks.mark_duplicates,
        args=(
            pypeliner.managed.TempInputFile('realigned_read_group_bam', 'read_group_id'),
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'tmp_dir' : pypeliner.managed.TempSpace('markdup_temp', 'read_group_id')
        }
    )

    return workflow

