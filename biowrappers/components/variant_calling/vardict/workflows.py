from pypeliner.workflow import Workflow

import pypeliner

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.variant_calling.utils as utils
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks


def create_vardict_single_sample_workflow(
        bam_file,
        ref_genome_fasta_file,
        indel_vcf_file,
        snv_vcf_file,
        chromosomes=default_chromosomes,
        min_allele_frequency=0.01,
        sample_name=None,
        split_size=int(1e7)):

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('config', 'regions'),
        value=utils.get_bam_regions(bam_file, split_size, chromosomes=chromosomes)
    )

    workflow.transform(
        name='run_vardict',
        axes=('regions',),
        ctx={'mem': 12, 'num_retry': 4, 'mem_retry_increment': 2},
        func=tasks.run_single_sample_vardict,
        args=(
            pypeliner.managed.InputFile(bam_file),
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.TempInputObj('config', 'regions'),
            pypeliner.managed.TempOutputFile('result.vcf', 'regions'),
        ),
        kwargs={
            'min_allele_frequency': min_allele_frequency,
            'sample_name': sample_name,
        },
    )

    workflow.transform(
        name='concatenate_vcf',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('result.vcf', 'regions'),
            pypeliner.managed.TempOutputFile('result.merged.vcf')
        )
    )

    workflow.transform(
        name='finalise_all_variants',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('result.merged.vcf'),
            pypeliner.managed.TempOutputFile('result.merged.vcf.gz')
        )
    )

    workflow.transform(
        name='filter_indels',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=tasks.filter_vcf,
        args=(
            pypeliner.managed.TempInputFile('result.merged.vcf.gz'),
            pypeliner.managed.TempOutputFile('indels.vcf'),
            'indel'
        )
    )

    workflow.transform(
        name='finalise_indels',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('indels.vcf'),
            pypeliner.managed.OutputFile(indel_vcf_file)
        )
    )

    workflow.transform(
        name='filter_snvs',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=tasks.filter_vcf,
        args=(
            pypeliner.managed.TempInputFile('result.merged.vcf.gz'),
            pypeliner.managed.TempOutputFile('snvs.vcf'),
            'snv'
        )
    )

    workflow.transform(
        name='finalise_snvs',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('snvs.vcf'),
            pypeliner.managed.OutputFile(snv_vcf_file)
        )
    )

    return workflow
