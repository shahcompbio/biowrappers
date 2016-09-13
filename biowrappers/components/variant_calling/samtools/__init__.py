import pypeliner

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.variant_calling.utils as utils
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks


def create_samtools_variant_calling_workflow(
        bam_file,
        ref_genome_fasta_file,
        indel_vcf_file,
        snv_vcf_file,
        chromosomes=default_chromosomes,
        split_size=int(1e7)
    ):
    
    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('regions_obj', 'regions'),
        value=utils.get_regions(bam_file, split_size, chromosomes=chromosomes)
    )
    
    workflow.transform(
        name='run_samtools_variant_calling',
        axes=('regions',),
        ctx={'mem' : 4},
        func=tasks.run_samtools_variant_calling,
        args=(
            pypeliner.managed.InputFile(bam_file),
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.TempOutputFile('variants.vcf', 'regions'),
        ),
        kwargs={
            'region': pypeliner.managed.TempInputObj('regions_obj', 'regions'),
        },
    )
    
    workflow.transform(
        name='merge_indels',
        ctx={'mem' : 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('variants.vcf', 'regions'),
            pypeliner.managed.TempOutputFile('indels.vcf'),
        ),
        kwargs={
            'variant_filter': 'indel',
        },
    )
    
    workflow.transform(
        name='merge_snvs',
        ctx={'mem' : 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('variants.vcf', 'regions'),
            pypeliner.managed.TempOutputFile('snvs.vcf'),
        ),
        kwargs={
            'variant_filter': 'snv',
        },
    )
    
    workflow.subworkflow(
        name='finalize_snvs',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempOutputFile('snvs.vcf'),
            pypeliner.managed.OutputFile(snv_vcf_file),
        )
    )
    
    workflow.subworkflow(
        name='finalize_indels',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempOutputFile('indels.vcf'),
            pypeliner.managed.OutputFile(snv_vcf_file),
        )
    )
    
    return workflow
    
    