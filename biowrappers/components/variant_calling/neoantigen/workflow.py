from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.variant_calling.utils as utils

import tasks


def create_hla_type_workflow(
    normal_bam_file,
    hla_type_file,
    config,
):
    workflow = Workflow()

    workflow.commandline(
        name='extract_chr6',
        args=(
            'samtools', 'view', '-bh', '-f', '2', '-F', '4',
            pypeliner.managed.InputFile(normal_bam_file),
            '6',
            '|',
            'samtools', 'collate', '-O', '-', pypeliner.managed.TempSpace('chr6_collate_temp'),
            '|',
            'samtools', 'bam2fq',
            '-1', pypeliner.managed.TempOutputFile('chr6_reads_1.fq'),
            '-2', pypeliner.managed.TempOutputFile('chr6_reads_2.fq'),
            '-',
        ),
    )

    workflow.transform(
        name='optitype',
        ctx={'mem': 24},
        func=tasks.run_optitype,
        args=(
            pypeliner.managed.TempInputFile('chr6_reads_1.fq'),
            pypeliner.managed.TempInputFile('chr6_reads_2.fq'),
            pypeliner.managed.OutputFile(hla_type_file),
            pypeliner.managed.TempSpace('optitype_temp'),
        )
    )

    return workflow
