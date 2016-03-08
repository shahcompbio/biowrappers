from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.variant_calling.utils as utils

import tasks

default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']

def create_mutect_workflow(
    normal_bam_file,
    tumour_bam_file,
    ref_genome_fasta_file,
    cosmic_vcf_file,
    dbsnp_vcf_file,
    out_file,
    chromosomes=default_chromosomes,
    split_size=int(1e7)):
    
    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('regions', 'regions'),
        value=utils.get_regions(tumour_bam_file, split_size, chromosomes=chromosomes)
    )
    
    workflow.transform(
        name='run_classify',
        axes=('regions',),
        ctx={'mem' : 6, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=tasks.run_mutect,
        args=(
            pypeliner.managed.InputFile(normal_bam_file),
            pypeliner.managed.InputFile(tumour_bam_file),
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.InputFile(cosmic_vcf_file),
            pypeliner.managed.InputFile(dbsnp_vcf_file),
            pypeliner.managed.TempInputObj('regions', 'regions'),
            pypeliner.managed.TempOutputFile('classified.vcf', 'regions')
        ),
    )

    workflow.transform(
        name='merge_vcf',
        ctx={'mem' : 8, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('classified.vcf', 'regions'),
            pypeliner.managed.TempOutputFile('merged.vcf'),
        )
    )
    
    workflow.transform(
        name='filter_snvs',
        ctx={'mem' : 2, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=vcf_tasks.filter_vcf,
        args=(
            pypeliner.managed.TempInputFile('merged.vcf'),
            pypeliner.managed.TempOutputFile('merged.filtered.vcf')
        )
    )
    
    workflow.subworkflow(
        name='finalise',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('merged.filtered.vcf'),
            pypeliner.managed.OutputFile(out_file)
        )
    )
    
    return workflow
    
