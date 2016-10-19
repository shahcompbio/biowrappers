from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.variant_calling.utils as utils
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks

sml_ctx = {'mem' : 2, 'num_retry' : 3, 'mem_retry_increment' : 2}
med_ctx = {'mem' : 4, 'num_retry' : 3, 'mem_retry_increment' : 4}

def create_snv_allele_counts_for_vcf_targets_workflow(
    bam_file,
    vcf_file,
    out_file,
    count_duplicates=False,
    min_bqual=0,
    min_mqual=0,
    split_size=int(1e7),
    table_name='snv_allele_counts'):
    
    workflow = Workflow()
    
    workflow.transform(
        name='get_regions',
        ret=pypeliner.managed.TempOutputObj('regions_obj', 'regions'),
        func=utils.get_vcf_regions,
        args=(
            pypeliner.managed.InputFile(vcf_file),
            split_size,
        ),
    )
    
    workflow.transform(
        name='get_snv_allele_counts_for_vcf_targets',
        axes=('regions',),
        ctx=med_ctx,
        func=tasks.get_snv_allele_counts_for_vcf_targets,
        args=(
            pypeliner.managed.InputFile(bam_file),
            pypeliner.managed.InputFile(vcf_file),
            pypeliner.managed.TempOutputFile('counts.h5', 'regions'),
            table_name
        ),
        kwargs={
            'count_duplicates' : count_duplicates,
            'min_bqual' : min_bqual,
            'min_mqual' : min_mqual,
            'region': pypeliner.managed.TempInputObj('regions_obj', 'regions'),
        }
    )
    
    workflow.transform(
        name='merge_snv_allele_counts',
        ctx=med_ctx,
        func=hdf5_tasks.concatenate_tables,
        args=(
            pypeliner.managed.TempInputFile('counts.h5', 'regions'),
            pypeliner.managed.OutputFile(out_file)
        ),
        kwargs={
            'in_memory' : False,
        }
    )
    
    return workflow

def create_snv_allele_counts_workflow(
    bam_file,
    out_file,
    table_name,
    chromosomes=utils.default_chromosomes,
    count_duplicates=False,
    min_bqual=0,
    min_mqual=0,
    report_non_variant_positions=True,
    report_zero_count_positions=False,
    split_size=int(1e7)):
    
    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('regions_obj', 'regions'),
        value=utils.get_bam_regions(bam_file, split_size, chromosomes=chromosomes)
    )
    
    workflow.transform(
        name='get_counts',
        axes=('regions',),
        ctx=med_ctx,
        func=tasks.get_snv_allele_counts_for_region,
        args=(
            pypeliner.managed.InputFile(bam_file),
            pypeliner.managed.TempOutputFile('counts.h5', 'regions'),
            pypeliner.managed.TempInputObj('regions_obj', 'regions'),
            table_name
        ),
        kwargs={
            'count_duplicates' : count_duplicates,
            'min_bqual' : min_bqual,
            'min_mqual' : min_mqual,
            'report_non_variant_positions' : report_non_variant_positions,
            'report_zero_count_positions' : report_zero_count_positions
        }
    )
    
    workflow.transform(
        name='concatenate_counts',
        ctx=med_ctx,
        func=hdf5_tasks.concatenate_tables,
        args=(
            pypeliner.managed.TempInputFile('counts.h5', 'regions'),
            pypeliner.managed.OutputFile(out_file)
        )
    )
    
    return workflow

def create_snv_variant_position_counts_workflow(
    normal_bam_file,
    tumour_bam_files,
    out_file,
    chromosomes=utils.default_chromosomes,
    count_duplicates=False,
    min_bqual=0,
    min_mqual=0,
    min_normal_depth=0,
    min_tumour_depth=0,
    min_variant_depth=0,
    report_strand_counts=False,
    split_size=int(1e7),
    table_group=''):
    
    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('regions_obj', 'regions'),
        value=utils.get_bam_regions(normal_bam_file, split_size, chromosomes=chromosomes)
    )
    
    tumour_input_files = {}
    
    for sample in tumour_bam_files:
        tumour_input_files[sample] = pypeliner.managed.InputFile(tumour_bam_files[sample])
    
    workflow.transform(
        name='get_counts',
        axes=('regions',),
        ctx=med_ctx,
        func=tasks.get_variant_position_counts,
        args=(
            pypeliner.managed.InputFile(normal_bam_file),
            tumour_input_files,
            pypeliner.managed.TempOutputFile('counts.h5', 'regions'),
            pypeliner.managed.TempInputObj('regions_obj', 'regions'),
            table_group
        ),
        kwargs={
            'count_duplicates' : count_duplicates,
            'min_bqual' : min_bqual,
            'min_mqual' : min_mqual,
            'min_normal_depth' : min_normal_depth,
            'min_tumour_depth' : min_tumour_depth,
            'min_variant_depth' : min_variant_depth,
            'report_strand_counts': report_strand_counts
        }
    )
    
    workflow.transform(
        name='concatenate_counts',
        axes=(),
        ctx=med_ctx,
        func=hdf5_tasks.concatenate_tables,
        args=(
            pypeliner.managed.TempInputFile('counts.h5', 'regions'),
            pypeliner.managed.OutputFile(out_file),
        )
    )
    
    return workflow
