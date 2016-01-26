from pypeliner.workflow import Workflow

import os
import pypeliner  

from biowrappers.components.utils import make_parent_directory
from biowrappers.components.variant_calling.utils import default_chromosomes  

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.variant_calling.annotated_db_status as annotated_db_status
import biowrappers.components.variant_calling.mappability as mappability
import biowrappers.components.variant_calling.museq as museq
import biowrappers.components.variant_calling.mutect as mutect
import biowrappers.components.variant_calling.snpeff as snpeff
import biowrappers.components.variant_calling.snv_allele_counts as snv_allele_counts
import biowrappers.components.variant_calling.strelka as strelka
import biowrappers.components.variant_calling.tri_nucleotide_context as tri_nucleotide_context
import biowrappers.components.variant_calling.vardict as vardict

def call_and_annotate_pipeline(
        config,
        normal_bam_path,
        tumour_bam_paths,
        ref_genome_fasta_file,
        raw_data_dir,
        results_file,
        chromosomes=default_chromosomes):
    
    workflow = Workflow()
    
    workflow.setobj(pypeliner.managed.OutputChunks('tumour_sample_id', axes_origin=[0, ]), tumour_bam_paths.keys())
    
    snv_vcf_files = {
        'museq' : pypeliner.managed.File(get_sample_out_file('museq', 'vcf.gz', raw_data_dir), 'tumour_sample_id'),
        'mutect' : pypeliner.managed.File(get_sample_out_file('mutect', 'vcf.gz', raw_data_dir), 'tumour_sample_id'),
        'strelka' : pypeliner.managed.File(get_sample_out_file('strelka', 'vcf.gz', raw_data_dir, sub_output='snv'), 'tumour_sample_id'),
        'vardict' : pypeliner.managed.File(get_sample_out_file('vardict', 'vcf.gz', raw_data_dir, sub_output='snv'), 'tumour_sample_id')
    }
    
    indel_vcf_files = {
        'strelka' : pypeliner.managed.File(get_sample_out_file('strelka', 'vcf.gz', raw_data_dir, sub_output='indel'), 'tumour_sample_id'),
        'vardict' : pypeliner.managed.File(get_sample_out_file('vardict', 'vcf.gz', raw_data_dir, sub_output='indel'), 'tumour_sample_id')
    }
    
    normal_bam_file = pypeliner.managed.File(normal_bam_path)
    
    tumour_bam_files = pypeliner.managed.File('tumour_bams', 'tumour_sample_id', fnames=tumour_bam_paths)
    
    ref_genome_fasta_file = pypeliner.managed.File(ref_genome_fasta_file)
    
    for prog in ('museq', 'mutect', 'strelka', 'vardict'):
        config[prog]['kwargs']['chromosomes'] = chromosomes

    workflow.subworkflow(
        name='museq',
        axes=('tumour_sample_id',),
        func=museq.museq_pipeline,
        args=(
            normal_bam_file.as_input(),
            tumour_bam_files.as_input(),
            ref_genome_fasta_file.as_input(),
            snv_vcf_files['museq'].as_output()
        ),
        kwargs=config['museq']['kwargs']
    )    
        
    workflow.subworkflow(
        name='mutect',
        axes=('tumour_sample_id',),
        func=mutect.mutect_pipeline,
        args=(
            normal_bam_file.as_input(),
            tumour_bam_files.as_input(),
            ref_genome_fasta_file.as_input(),
            config['annotate_cosmic_status']['db_vcf_file'],
            config['annotate_dbsnp_status']['db_vcf_file'],
            snv_vcf_files['mutect'].as_output()
        ),
        kwargs=config['mutect']['kwargs']
    )

    workflow.subworkflow(
        name='strelka',
        axes=('tumour_sample_id',),
        func=strelka.strelka_pipeline,
        args=(
            normal_bam_file.as_input(),
            tumour_bam_files.as_input(),
            ref_genome_fasta_file.as_input(),
            indel_vcf_files['strelka'].as_output(),
            snv_vcf_files['strelka'].as_output()
        ),
        kwargs=config['strelka']['kwargs']
    )
      
    workflow.subworkflow(
        name='vardict',
        axes=('tumour_sample_id',),
        func=vardict.vardict_pipeline,
        args=(
            normal_bam_file.as_input(),
            tumour_bam_files.as_input(),
            ref_genome_fasta_file.as_input(),
            indel_vcf_files['vardict'].as_output(),
            snv_vcf_files['vardict'].as_output()
        ),
        kwargs=config['vardict']['kwargs']
    )
     
    workflow.transform(
        name='merge_indels',
        ctx={'mem' : 8},
        func=vcf_tasks.merge_vcfs,
        args=(
            [x.as_input() for x in indel_vcf_files.values()],
            pypeliner.managed.TempOutputFile('all.indel.vcf')
        )
    )
      
    workflow.subworkflow(
        name='finalise_indels',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('all.indel.vcf'),
            pypeliner.managed.TempOutputFile('all.indel.vcf.gz')
        )
    )
      
    workflow.transform(
        name='merge_snvs',
        ctx={'mem' : 8},
        func=vcf_tasks.merge_vcfs,
        args=(
            [x.as_input() for x in snv_vcf_files.values()],
            pypeliner.managed.TempOutputFile('all.snv.vcf')
        )
    )
      
    workflow.subworkflow(
        name='finalise_snvs',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('all.snv.vcf'),
            pypeliner.managed.TempOutputFile('all.snv.vcf.gz')
        )
    )
      
    workflow.subworkflow(
        name='snpeff_snvs',
        func=snpeff.snpeff_pipeline,
        args=(
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('snpeff.h5'),
        ),
        kwargs=config['snpeff']['kwargs']
    )

    workflow.subworkflow(
        name='annotate_snv_cosmic_status',
        func=annotated_db_status.vcf_db_annotation_pipeline,
        args=(
            pypeliner.managed.InputFile(config['annotate_cosmic_status']['db_vcf_file']),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('cosmic.h5'),
        ),
        kwargs=config['annotate_cosmic_status']['kwargs']
    )
    
    workflow.subworkflow(
        name='annotate_snv_dbsnp_status',
        func=annotated_db_status.vcf_db_annotation_pipeline,
        args=(
            pypeliner.managed.InputFile(config['annotate_dbsnp_status']['db_vcf_file']),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('dbsnp.h5'),
        ),
        kwargs=config['annotate_dbsnp_status']['kwargs']
    )    
    
    workflow.subworkflow(
        name='snv_mappability',
        func=mappability.vcf_mappability_annotation_pipeline,
        args=(
            pypeliner.managed.InputFile(config['mappability']['mappability_file']),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('mappability.h5')
        ),
        kwargs=config['mappability']['kwargs']
    )
     
    workflow.subworkflow(
        name='snv_tri_nucleotide_context',
        func=tri_nucleotide_context.vcf_tric_nucleotide_annotation_pipeline,
        args=(
            ref_genome_fasta_file.as_input(),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('tri_nucleotide_context.h5')
        ),
        kwargs=config['tri_nucleotide_context']['kwargs']
    )
    
    normal_snv_counts_kwargs = config['snv_counts']['kwargs'].copy()
    normal_snv_counts_kwargs['table_name'] = '/snv/counts/normal'
    
    workflow.subworkflow(
        name='normal_snv_counts',
        func=snv_allele_counts.snv_allele_counts_for_vcf_targets_pipeline,
        args=(
            normal_bam_file.as_input(),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('normal_counts.h5'),
        ),
        kwargs=normal_snv_counts_kwargs
    )
  
    tumour_snv_counts_kwargs = config['snv_counts']['kwargs'].copy()
    tumour_snv_counts_kwargs['table_name'] = pypeliner.managed.Template('/snv/counts/{tumour_sample_id}', 'tumour_sample_id')
    
    workflow.subworkflow(
        name='tumour_snv_counts',
        axes=('tumour_sample_id',),
        func=snv_allele_counts.snv_allele_counts_for_vcf_targets_pipeline,
        args=(
            tumour_bam_files.as_input(),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('tumour_counts.h5', 'tumour_sample_id')
        ),
        kwargs=tumour_snv_counts_kwargs
    )
    
    workflow.transform(
        name='build_results_file',
        ctx={'mem' : 2},
        func=hdf5_tasks.concatenate_tables,
        args=((
               pypeliner.managed.TempInputFile('cosmic.h5'),
               pypeliner.managed.TempInputFile('snpeff.h5'),
               pypeliner.managed.TempInputFile('mappability.h5'),
               pypeliner.managed.TempInputFile('tri_nucleotide_context.h5'),
               pypeliner.managed.TempInputFile('normal_counts.h5'),
               pypeliner.managed.TempInputFile('tumour_counts.h5', 'tumour_sample_id')
            ),
            pypeliner.managed.OutputFile(results_file)
        ),
        
    )
    
    return workflow

def get_sample_out_file(cmd, ext, out_dir, sub_output=None):
    if sub_output is not None:
        cmd = os.path.join(cmd, sub_output)
        
    out_file = os.path.join(out_dir, cmd, '{{tumour_sample_id}}.{0}'.format(ext))
    
    make_parent_directory(out_file)
    
    return out_file
