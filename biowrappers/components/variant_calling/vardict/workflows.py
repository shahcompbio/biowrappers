from pypeliner.workflow import Workflow

import pypeliner

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.variant_calling.utils as utils
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks


def create_vardict_single_sample_workflow(
        bam_file,
        ref_genome_fasta_file,
        out_file,
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
        name='compress_tmp',
        axes=('regions',),
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=vcf_tasks.compress_vcf,
        args=(
            pypeliner.managed.TempInputFile('result.vcf', 'regions'),
            pypeliner.managed.TempOutputFile('result.vcf.gz', 'regions'),
        ),
        kwargs={
            'index_file': pypeliner.managed.TempOutputFile('result.vcf.gz.tbi', 'regions'),
        }
    )
    workflow.transform(
        name='concatenate_vcf',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('result.vcf.gz', 'regions'),
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'bcf_index_file': pypeliner.managed.OutputFile(out_file + '.csi'),
            'vcf_index_file': pypeliner.managed.OutputFile(out_file + '.tbi'),
        },
    )
    return workflow
