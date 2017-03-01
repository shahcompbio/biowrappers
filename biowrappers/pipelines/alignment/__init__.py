import pypeliner

import biowrappers.components.io.bam.tasks
import biowrappers.components.io.fastq.tasks
import biowrappers.components.alignment.bwa.tasks
import biowrappers.pipelines.realignment.tasks


def alignment_pipeline(
        config,
        fastq_1,
        fastq_2,
        out_file):

    ref_genome = pypeliner.managed.InputFile(config['ref_genome']['file'])

    read_group_config = config.get('read_group', {})

    if 'ID' not in read_group_config:
        read_group_config['ID'] = hash(fastq_1 + fastq_2) % int(1e6)

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('read_group_config'),
        value=read_group_config,
    )

    workflow.transform(
        name='split_fastq_1',
        ctx={'mem': 4},
        func=biowrappers.components.io.fastq.tasks.split_fastq,
        args=(
            pypeliner.managed.InputFile(fastq_1),
            pypeliner.managed.TempOutputFile('read_1', 'split'),
            config['split_size'],
        ),
    )

    workflow.transform(
        name='split_fastq_2',
        ctx={'mem': 4},
        func=biowrappers.components.io.fastq.tasks.split_fastq,
        args=(
            pypeliner.managed.InputFile(fastq_2),
            pypeliner.managed.TempOutputFile('read_2', 'split', axes_origin=[]),
            config['split_size'],
        ),
    )

    workflow.transform(
        name='aln_read_1',
        axes=('split',),
        ctx={'mem': 6},
        func=biowrappers.components.alignment.bwa.tasks.run_aln,
        args=(
            pypeliner.managed.TempInputFile('read_1', 'split'),
            ref_genome,
            pypeliner.managed.TempOutputFile('read_1.sai', 'split'),
        ),
    )

    workflow.transform(
        name='aln_read_2',
        axes=('split',),
        ctx={'mem': 6},
        func=biowrappers.components.alignment.bwa.tasks.run_aln,
        args=(
            pypeliner.managed.TempInputFile('read_2', 'split'),
            ref_genome,
            pypeliner.managed.TempOutputFile('read_2.sai', 'split'),
        ),
    )

    workflow.transform(
        name='sampe',
        axes=('split',),
        ctx={'mem': 6},
        func=biowrappers.components.alignment.bwa.tasks.run_sampe,
        args=(
            pypeliner.managed.TempInputFile('read_1', 'split'),
            pypeliner.managed.TempInputFile('read_2', 'split'),
            pypeliner.managed.TempInputFile('read_1.sai', 'split'),
            pypeliner.managed.TempInputFile('read_2.sai', 'split'),
            ref_genome,
            pypeliner.managed.TempOutputFile('aligned.bam', 'split'),
        ),
        kwargs={
            'read_group_info': pypeliner.managed.TempInputObj('read_group_config'),
        },
    )

    workflow.transform(
        name='sort',
        axes=('split',),
        ctx={'mem': 4},
        func=biowrappers.components.io.bam.tasks.sort,
        args=(
            pypeliner.managed.TempInputFile('aligned.bam', 'split'),
            pypeliner.managed.TempOutputFile('sorted.bam', 'split'),
        ),
    )

    workflow.transform(
        name='write_header_file',
        ctx={'local': True},
        func=biowrappers.pipelines.realignment.tasks.write_header_file,
        args=(
            pypeliner.managed.TempInputFile('sorted.bam', 'split'),
            pypeliner.managed.TempOutputFile('header.sam'),
            config['ref_genome']['header']
        ),
    )

    workflow.transform(
        name='merge',
        ctx={'mem': 4},
        func=biowrappers.components.io.bam.tasks.merge,
        args=(
            pypeliner.managed.TempInputFile('sorted.bam', 'split'),
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'header_file': pypeliner.managed.TempInputFile('header.sam'),
        },
    )

    return workflow
