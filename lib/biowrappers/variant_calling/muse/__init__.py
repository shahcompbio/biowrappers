from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.io.vcf.tasks as vcf_tasks
import biowrappers.variant_calling.utils as utils

import tasks

default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']

def muse_pipeline(
    normal_bam_file,
    tumour_bam_file,
    ref_genome_fasta_file,
    out_file,
    chromosomes=default_chromosomes,
    exome=False,
    split_size=int(1e6)):
    
    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('muse_config', 'regions'),
        value=utils.get_regions(tumour_bam_file, split_size, chromosomes=chromosomes)
    )
    
    workflow.transform(
        name='run_muse_call',
        axes=('regions',),
        ctx={'mem' : 4},
        func=tasks.run_muse_call,
        args=(
            pypeliner.managed.InputFile(normal_bam_file),
            pypeliner.managed.InputFile(tumour_bam_file),
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.TempOutputFile('muse.split.MuSE.txt', 'regions')
        ),
        kwargs={'region' : pypeliner.managed.TempInputObj('muse_config', 'regions')}
    )
    
    workflow.transform(
        name='run_muse_sump',
        axes=('regions',),
        ctx={'mem' : 4},
        func=tasks.run_muse_sump,
        args=(
            pypeliner.managed.TempInputFile('muse.split.MuSE.txt', 'regions'),
            pypeliner.managed.TempOutputFile('muse.split.MuSE.sump.vcf', 'regions')
        ),
        kwargs={'exome' : exome}
    )
      
    workflow.transform(
        name='merge_muse_vcf',
        ctx={'mem' : 8},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('muse.split.MuSE.sump.vcf', 'regions'),
            pypeliner.managed.TempOutputFile('muse.merged.vcf')
        )
    )
    
    workflow.transform(
        name='filter_snvs',
        ctx={'mem' : 2},
        func=vcf_tasks.filter_vcf,
        args=(
            pypeliner.managed.TempInputFile('muse.merged.vcf'),
            pypeliner.managed.TempOutputFile('muse.merged.filtered.vcf')
        )
    )
    
    workflow.subworkflow(
        name='finalise',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('muse.merged.filtered.vcf'),
            pypeliner.managed.OutputFile(out_file)
        )
    )
    
    return workflow
    