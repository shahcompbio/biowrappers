from pypeliner.workflow import Workflow

import os
import pypeliner  

from biowrappers.components.utils import make_parent_directory
from biowrappers.components.variant_calling.utils import default_chromosomes  

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.variant_calling.annotated_db_status as annotated_db_status
import biowrappers.components.variant_calling.mappability as mappability
import biowrappers.components.variant_calling.nuseq as nuseq
import biowrappers.components.variant_calling.mutect as mutect
import biowrappers.components.variant_calling.snpeff as snpeff
import biowrappers.components.variant_calling.snv_allele_counts as snv_allele_counts
import biowrappers.components.variant_calling.strelka as strelka
import biowrappers.components.variant_calling.tri_nucleotide_context as tri_nucleotide_context
import biowrappers.components.variant_calling.vardict as vardict
    
default_ctx = {'mem' : 4, 'num_retry' : 3, 'mem_retry_increment' : 2}

big_mem_ctx = {'mem' : 8, 'num_retry' : 3, 'mem_retry_increment' : 2}

vcf_score_callbacks = {
    'indel' : {
        'strelka' : lambda record: record.INFO['QSI'],
        'vardict' : None
    },
    'snv': {
        'mutect' : None,
        'nuseq' : lambda record: record.INFO['PS'],
        'strelka' : lambda record: record.INFO['QSI'],
        'vardict' : None
    }
}

def call_and_annotate_pipeline(
        config,
        normal_bam_path,
        tumour_bam_paths,
        raw_data_dir,
        results_file,
        chromosomes=default_chromosomes):

    workflow = Workflow()
    
    workflow.setobj(pypeliner.managed.OutputChunks('tumour_sample_id', axes_origin=[0, ]), tumour_bam_paths.keys())
    
    indel_vcf_files = {}
    
    snv_vcf_files = {}
    
    indel_progs = ('strelka', 'vardict')
    
    indel_progs = [x for x in indel_progs if x in config]
    
    snv_progs = ('nuseq', 'mutect', 'strelka', 'vardict')
    
    snv_progs = [x for x in snv_progs if x in config]
    
    for prog in snv_progs:
        config[prog]['kwargs']['chromosomes'] = chromosomes
        
        snv_vcf_files[prog] = pypeliner.managed.File(
            get_sample_out_file(prog, 'vcf.gz', raw_data_dir),
            'tumour_sample_id'
        )
        
        if prog in indel_progs:
            indel_vcf_files[prog] = pypeliner.managed.File(
                get_sample_out_file(prog, 'vcf.gz', raw_data_dir, sub_output='indel'),
                'tumour_sample_id'
            )
    
    if 'nuseq_multi_sample' in config:
        config['nuseq_multi_sample']['kwargs']['chromosomes'] = chromosomes
        
        snv_vcf_files['nuseq_multi_sample'] = pypeliner.managed.File(
            get_sample_out_file('nuseq_multi_sample', 'vcf.gz', raw_data_dir).format(tumour_sample_id='all')
        )
    
    normal_bam_file = pypeliner.managed.File(normal_bam_path)
    
    tumour_bam_files = pypeliner.managed.File('tumour_bams', 'tumour_sample_id', fnames=tumour_bam_paths)
    
    ref_genome_fasta_file = pypeliner.managed.File(config['databases']['ref_genome']['local_path'])
    
    #===================================================================================================================
    # Multi sample calling
    #===================================================================================================================
    if 'nuseq_multi_sample' in config:
        workflow.subworkflow(
            name='nuseq_multi_sample',
            axes=(),
            func=nuseq.create_nuseq_classify_workflow,
            args=(
                normal_bam_file.as_input(),
                [pypeliner.managed.InputFile(x) for x in tumour_bam_paths.values()],
                ref_genome_fasta_file.as_input(),
                snv_vcf_files['nuseq_multi_sample'].as_output()
            ),
            kwargs=config['nuseq_multi_sample']['kwargs']
        )
            
        workflow.transform(
            name='convert_nuseq_multi_sample_vcf_to_hdf5',
            axes=(),
            ctx=default_ctx,
            func=vcf_tasks.convert_vcf_to_hdf5,
            args=(
                snv_vcf_files['nuseq_multi_sample'].as_input(),
                pypeliner.managed.TempOutputFile('nuseq_multi_sample.h5'),
                '/snv/vcf/nuseq_multi_sample/all',
            ),
            kwargs={
                'score_callback' : vcf_score_callbacks['nuseq']
            }
        )
    #===================================================================================================================
    # Single sample calling
    #===================================================================================================================
    if 'nuseq' in config:
        workflow.subworkflow(
            name='nuseq',
            axes=('tumour_sample_id',),
            func=nuseq.create_nuseq_classify_workflow,
            args=(
                normal_bam_file.as_input(),
                tumour_bam_files.as_input(),
                ref_genome_fasta_file.as_input(),
                snv_vcf_files['nuseq'].as_output()
            ),
            kwargs=config['nuseq']['kwargs']
        )
                
    if 'mutect' in config:
        workflow.subworkflow(
            name='mutect',
            axes=('tumour_sample_id',),
            func=mutect.create_mutect_workflow,
            args=(
                normal_bam_file.as_input(),
                tumour_bam_files.as_input(),
                ref_genome_fasta_file.as_input(),
                config['databases']['cosmic']['local_path'],
                config['databases']['dbsnp']['local_path'],
                snv_vcf_files['mutect'].as_output()
            ),
            kwargs=config['mutect']['kwargs']
        )
    
    if 'strelka' in config:
        workflow.subworkflow(
            name='strelka',
            axes=('tumour_sample_id',),
            func=strelka.create_strelka_workflow,
            args=(
                normal_bam_file.as_input(),
                tumour_bam_files.as_input(),
                ref_genome_fasta_file.as_input(),
                indel_vcf_files['strelka'].as_output(),
                snv_vcf_files['strelka'].as_output()
            ),
            kwargs=config['strelka']['kwargs']
        )
    
    if 'vardict' in config:
        workflow.subworkflow(
            name='vardict',
            axes=('tumour_sample_id',),
            func=vardict.create_vardict_workflow,
            args=(
                normal_bam_file.as_input(),
                tumour_bam_files.as_input(),
                ref_genome_fasta_file.as_input(),
                indel_vcf_files['vardict'].as_output(),
                snv_vcf_files['vardict'].as_output()
            ),
            kwargs=config['vardict']['kwargs']
        )
    
    #===================================================================================================================
    # Indel annotation
    #===================================================================================================================
    for prog in indel_progs:
        workflow.transform(
            name='convert_{0}_indel_vcf_to_hdf5'.format(prog),
            axes=('tumour_sample_id',),
            ctx=default_ctx,
            func=vcf_tasks.convert_vcf_to_hdf5,
            args=(
                indel_vcf_files[prog].as_input(),
                pypeliner.managed.TempOutputFile('indel_{0}.h5'.format(prog), 'tumour_sample_id'),
                pypeliner.managed.Template('/indel/vcf/{prog}/{{tumour_sample_id}}'.format(prog=prog), 'tumour_sample_id')
            ),
            kwargs={
                'score_callback' : vcf_score_callbacks['indel'][prog] 
            }
        )
    
    workflow.transform(
        name='merge_indels',
        ctx=big_mem_ctx,
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
    
    #===================================================================================================================
    # SNV
    #===================================================================================================================
    for prog in snv_progs:
        workflow.transform(
            name='convert_{0}_snv_vcf_to_hdf5'.format(prog),
            axes=('tumour_sample_id',),
            ctx=default_ctx,
            func=vcf_tasks.convert_vcf_to_hdf5,
            args=(
                snv_vcf_files[prog].as_input(),
                pypeliner.managed.TempOutputFile('snv_{0}.h5'.format(prog), 'tumour_sample_id'),
                pypeliner.managed.Template('/snv/vcf/{prog}/{{tumour_sample_id}}'.format(prog=prog), 'tumour_sample_id')
            ),
            kwargs={
                'score_callback' : vcf_score_callbacks['snv'][prog] 
            }
        )
        
    workflow.transform(
        name='merge_snvs',
        ctx=big_mem_ctx,
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
        name='annotate_snvs', 
        axes=(), 
        func=create_annotation_workflow, 
        args=(
            config,
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('snv_annotations.h5'),
        ),
        kwargs={
            'variant_type' : 'snv'
        }
    )

    workflow.subworkflow(
        name='normal_snv_counts',
        func=snv_allele_counts.create_snv_allele_counts_for_vcf_targets_workflow,
        args=(
            normal_bam_file.as_input(),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('normal_counts.h5'),
        ),
        kwargs=get_kwargs(config['snv_counts']['kwargs'], '/snv/counts/normal')
    )

    workflow.subworkflow(
        name='tumour_snv_counts',
        axes=('tumour_sample_id',),
        func=snv_allele_counts.create_snv_allele_counts_for_vcf_targets_workflow,
        args=(
            tumour_bam_files.as_input(),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.TempOutputFile('tumour_counts.h5', 'tumour_sample_id')
        ),
        kwargs=get_kwargs(config['snv_counts']['kwargs'], pypeliner.managed.Template('/snv/counts/{tumour_sample_id}', 'tumour_sample_id'))
    )
    
    tables = [
        pypeliner.managed.TempInputFile('snv_annotations.h5'),
        pypeliner.managed.TempInputFile('normal_counts.h5'),
        pypeliner.managed.TempInputFile('tumour_counts.h5', 'tumour_sample_id'),
    ]

    for prog in indel_progs:
        tables.append(pypeliner.managed.TempInputFile('indel_{0}.h5'.format(prog), 'tumour_sample_id'))
            
    for prog in snv_progs:
        tables.append(pypeliner.managed.TempInputFile('snv_{0}.h5'.format(prog), 'tumour_sample_id'))
    
    if 'nuseq_multi_sample' in config:
        tables.append(pypeliner.managed.TempInputFile('nuseq_multi_sample.h5'))
    
    workflow.transform(
        name='build_results_file',
        ctx=default_ctx,
        func=hdf5_tasks.concatenate_tables,
        args=(
            tables,
            pypeliner.managed.OutputFile(results_file)
        ),
        
    )
    
    return workflow

def create_annotation_workflow(config, in_vcf_file, out_file, variant_type='snv'):
    
    annotators = (
        'cosmic_status', 
        'dbsnp_status',
        'mappability',
        'snpeff',
        'tri_nucleotide_context'
    )
    
    kwargs = {}
    
    temp_result_files = {}
    
    for a in annotators:
        kwargs[a] = get_kwargs(config[a]['kwargs'], '/{0}/{1}'.format(variant_type, a))
        
        temp_result_files[a] = pypeliner.managed.TempFile('{0}.h5'.format(a))
    
    workflow= Workflow()
    
    workflow.subworkflow(
        name='snpeff',
        func=snpeff.create_snpeff_annotation_workflow,
        args=(
            config['databases']['snpeff']['db'],
            pypeliner.managed.InputFile(in_vcf_file),
            temp_result_files['snpeff'].as_output(),
        ),
        kwargs=kwargs['snpeff']
    )

    workflow.subworkflow(
        name='cosmic_status',
        func=annotated_db_status.create_vcf_db_annotation_workflow,
        args=(
            pypeliner.managed.InputFile(config['databases']['cosmic']['local_path']),
            pypeliner.managed.InputFile(in_vcf_file),
            temp_result_files['cosmic_status'].as_output(),
        ),
        kwargs=kwargs['cosmic_status']
    )
    
    workflow.subworkflow(
        name='dbsnp_status',
        func=annotated_db_status.create_vcf_db_annotation_workflow,
        args=(
            pypeliner.managed.InputFile(config['databases']['dbsnp']['local_path']),
            pypeliner.managed.InputFile(in_vcf_file),
            temp_result_files['dbsnp_status'].as_output(),
        ),
        kwargs=kwargs['dbsnp_status']
    )    
    
    workflow.subworkflow(
        name='mappability',
        func=mappability.create_vcf_mappability_annotation_workflow,
        args=(
            pypeliner.managed.InputFile(config['databases']['mappability']['local_path']),
            pypeliner.managed.InputFile(in_vcf_file),
            temp_result_files['mappability'].as_output(),
        ),
        kwargs=kwargs['mappability']
    )
     
    workflow.subworkflow(
        name='snv_tri_nucleotide_context',
        func=tri_nucleotide_context.create_vcf_tric_nucleotide_annotation_workflow,
        args=(
            pypeliner.managed.InputFile(config['databases']['ref_genome']['local_path']),
            pypeliner.managed.InputFile(in_vcf_file),
            temp_result_files['tri_nucleotide_context'].as_output(),
        ),
        kwargs=kwargs['snv_tri_nucleotide_context']
    )
    
    workflow.transform(
        name='build_results_file',
        ctx=default_ctx,
        func=hdf5_tasks.concatenate_tables,
        args=(
            [x.as_input() for x in temp_result_files.values()],
            pypeliner.managed.OutputFile(out_file)
        ),
    )
    
    return workflow
    
def get_sample_out_file(cmd, ext, out_dir, sub_output=None):
    if sub_output is not None:
        cmd = os.path.join(cmd, sub_output)
        
    out_file = os.path.join(out_dir, cmd, '{{tumour_sample_id}}.{0}'.format(ext))
    
    make_parent_directory(out_file)
    
    return out_file

def get_kwargs(config, table_name):
    
    config = config.copy()
    
    config['table_name'] = table_name
    
    return config
    