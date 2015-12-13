from collections import OrderedDict

import os
import pypeliner
import pysam

import biowrappers.variant_calling.utils as utils
import biowrappers.io.vcf.tasks as vcf_tasks
import tasks

default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']

def mutation_seq_pipeline(
    sch,
    normal_bam_file,
    tumour_bam_file,
    ref_genome_fasta_file,
    out_file,
    install_dir,
    chromosomes=default_chromosomes,
    somatic_threshold=0.85,
    split_size=int(1e6)):
    
    sch.setobj(
        pypeliner.managed.TempOutputObj('museq_config', 'regions'),
        utils.get_regions(tumour_bam_file, split_size, chromosomes=chromosomes)
    )
    
    sch.transform(
        'run_museq',
        ('regions',),
        {'mem' : 4},
        tasks.run_museq,
        None,
        install_dir,
        pypeliner.managed.InputFile(normal_bam_file),
        pypeliner.managed.InputFile(tumour_bam_file),
        pypeliner.managed.InputFile(ref_genome_fasta_file),
        pypeliner.managed.TempOutputFile('museq.split', 'regions'),
        region=pypeliner.managed.TempInputObj('museq_config', 'regions')
    )
    
    sch.transform(
        'filter_museq_calls',
        ('regions',),
        {'mem' : 4},
        tasks.filter_snvs_by_threshold,
        None,
        pypeliner.managed.TempInputFile('museq.split', 'regions'),
        pypeliner.managed.TempOutputFile('museq.split.filtered', 'regions'),
        threshold=somatic_threshold
    )
    
    sch.transform(
        'merge_museq_vcf',
        (),
        {'mem' : 8},
        vcf_tasks.concatenate_vcf,
        None,
        pypeliner.managed.TempInputFile('museq.split.filtered', 'regions'),
        pypeliner.managed.TempOutputFile('museq.merged')
    )
    
    sch.transform(
        'compress_vcf',
        (),
        {'mem' : 2},
        vcf_tasks.compress_vcf,
        None,
        pypeliner.managed.TempInputFile('museq.merged'),
        pypeliner.managed.OutputFile(out_file)
    )
    
    sch.transform(
        'index_vcf',
        (),
        {'mem' : 2},
        vcf_tasks.index_vcf,
        None,
        pypeliner.managed.InputFile(out_file),
        pypeliner.managed.OutputFile(out_file + '.tbi')
    )
    