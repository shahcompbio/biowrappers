import pypeliner

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.variant_calling.utils as utils
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks

def samtools_variant_calling_pipeline(
        sch,
        bam_file,
        ref_genome_fasta_file,
        indel_vcf_file,
        snv_vcf_file,
        chromosomes=default_chromosomes,
        split_size=int(1e6)
    ):
    
    sch.setobj(
        pypeliner.managed.TempOutputObj('regions_obj', 'regions'),
        utils.get_regions(bam_file, split_size, chromosomes=chromosomes)
    )
    
    sch.transform(
        'run_samtools_variant_calling',
        ('regions',),
        {'mem' : 4},
        tasks.run_samtools_variant_calling,
        None,
        pypeliner.managed.InputFile(bam_file),
        pypeliner.managed.InputFile(ref_genome_fasta_file),
        pypeliner.managed.TempOutputFile('variants.vcf', 'regions'),
        region=pypeliner.managed.TempInputObj('regions_obj', 'regions')
    )
    
    sch.transform(
        'merge_indels',
        (),
        {'mem' : 2},
        vcf_tasks.concatenate_vcf,
        None,
        pypeliner.managed.TempInputFile('variants.vcf', 'regions'),
        pypeliner.managed.TempOutputFile('indels.vcf'),
        variant_filter='indel'
    )
    
    sch.transform(
        'merge_snvs',
        (),
        {'mem' : 2},
        vcf_tasks.concatenate_vcf,
        None,
        pypeliner.managed.TempInputFile('variants.vcf', 'regions'),
        pypeliner.managed.TempOutputFile('snvs.vcf'),
        variant_filter='snv'
    )
    
    vcf_tasks.finalise_vcf(sch, 'indel', 'indels.vcf', indel_vcf_file)
    
    vcf_tasks.finalise_vcf(sch, 'snv', 'snvs.vcf', snv_vcf_file)
    
    