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
    
    workflow = pypeliner.workflow.Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('regions_obj', 'regions'),
        value=utils.get_bam_regions(bam_file, split_size, chromosomes=chromosomes)
    )
    
    workflow.transform(
        name='run_samtools_variant_calling',
        axes=('regions',),
        ctx={'mem' : 4},
        func=tasks.run_samtools_variant_calling,
        args=(
            pypeliner.managed.InputFile(bam_file),
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.TempOutputFile('variants.vcf.gz', 'regions'),
        ),
        kwargs={
            'region': pypeliner.managed.TempInputObj('regions_obj', 'regions'),
        },
    )
    
    workflow.transform(
        name='concatenate_variants',
        ctx={'mem' : 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('variants.vcf.gz', 'regions'),
            pypeliner.managed.TempOutputFile('variants.vcf.gz'),
        ),
    )
    
    workflow.transform(
        name='extract_indels',
        ctx={'mem' : 2},
        func=vcf_tasks.extract_variant_type,
        args=(
            pypeliner.managed.TempInputFile('variants.vcf.gz'),
            pypeliner.managed.OutputFile(indel_vcf_file),
            'indels',
        ),
    )
    
    workflow.transform(
        name='extract_snvs',
        ctx={'mem' : 2},
        func=vcf_tasks.extract_variant_type,
        args=(
            pypeliner.managed.TempInputFile('variants.vcf.gz'),
            pypeliner.managed.OutputFile(snv_vcf_file),
            'snps',
        ),
    )
    
    return workflow
    
    
