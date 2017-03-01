import pypeliner

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.variant_calling.utils as utils
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks


def create_samtools_variant_calling_workflow(
    bam_file,
    ref_genome_fasta_file,
    vcf_file,
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
        ctx={'mem': 4},
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
        ctx={'mem': 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('variants.vcf.gz', 'regions'),
            pypeliner.managed.OutputFile(vcf_file),
        ),
    )

    return workflow
