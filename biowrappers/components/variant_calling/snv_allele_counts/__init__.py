import pypeliner
import pypeliner.managed as mgd
import biowrappers

default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


sml_ctx = {'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2}
med_ctx = {'mem': 4, 'num_retry': 3, 'mem_retry_increment': 4}


def create_snv_allele_counts_for_vcf_targets_workflow(
        bam_file,
        vcf_file,
        out_file,
        chromosomes=default_chromosomes,
        count_duplicates=False,
        hdf5_output=True,
        min_bqual=0,
        min_mqual=0,
        split_size=int(1e7),
        table_name='snv_allele_counts',
        vcf_to_bam_chrom_map=None):

    if hdf5_output:
        merged_file = mgd.File(out_file)

    else:
        merged_file = mgd.TempFile('merged.h5')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='get_regions',
        ret=mgd.TempOutputObj('regions_obj', 'regions'),
        func='biowrappers.components.variant_calling.utils.get_vcf_regions',
        args=(
            mgd.InputFile(vcf_file),
            split_size,
        ),
        kwargs={
            'chromosomes': chromosomes,
        },
    )

    workflow.transform(
        name='get_snv_allele_counts_for_vcf_targets',
        axes=('regions',),
        ctx=med_ctx,
        func='biowrappers.components.snv_allele_counts.tasks.get_snv_allele_counts_for_vcf_targets',
        args=(
            mgd.InputFile(bam_file),
            mgd.InputFile(vcf_file),
            mgd.TempOutputFile('counts.h5', 'regions'),
            table_name
        ),
        kwargs={
            'count_duplicates': count_duplicates,
            'min_bqual': min_bqual,
            'min_mqual': min_mqual,
            'region': mgd.TempInputObj('regions_obj', 'regions'),
            'vcf_to_bam_chrom_map': vcf_to_bam_chrom_map,
        }
    )

    workflow.transform(
        name='merge_snv_allele_counts',
        ctx=med_ctx,
        func='biowrappers.components.io.hdf5.tasks.concatenate_tables',
        args=(
            mgd.TempInputFile('counts.h5', 'regions'),
            merged_file.as_output(),
        ),
        kwargs={
            'in_memory': False,
        }
    )

    if not hdf5_output:
        workflow.transform(
            name='convert_to_tsv',
            ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
            func='biowrappers.components.io.hdf5.tasks.convert_hdf5_to_tsv',
            args=(
                merged_file.as_input(),
                table_name,
                mgd.OutputFile(out_file),
            ),
            kwargs={
                'compress': True,
            }
        )

    return workflow


def create_snv_allele_counts_workflow(
        bam_file,
        out_file,
        table_name,
        chromosomes=default_chromosomes,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        report_non_variant_positions=True,
        report_zero_count_positions=False,
        split_size=int(1e7)):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.TempOutputObj('regions_obj', 'regions'),
        value=biowrappers.components.variant_calling.utils.get_bam_regions(bam_file, split_size, chromosomes=chromosomes)
    )

    workflow.transform(
        name='get_counts',
        axes=('regions',),
        ctx=med_ctx,
        func='biowrappers.components.snv_allele_counts.tasks.get_snv_allele_counts_for_region',
        args=(
            mgd.InputFile(bam_file),
            mgd.TempOutputFile('counts.h5', 'regions'),
            mgd.TempInputObj('regions_obj', 'regions'),
            table_name
        ),
        kwargs={
            'count_duplicates': count_duplicates,
            'min_bqual': min_bqual,
            'min_mqual': min_mqual,
            'report_non_variant_positions': report_non_variant_positions,
            'report_zero_count_positions': report_zero_count_positions
        }
    )

    workflow.transform(
        name='concatenate_counts',
        ctx=med_ctx,
        func='biowrappers.components.io.hdf5.tasks.concatenate_tables',
        args=(
            mgd.TempInputFile('counts.h5', 'regions'),
            mgd.OutputFile(out_file)
        )
    )

    return workflow


def create_snv_variant_position_counts_workflow(
        normal_bam_file,
        tumour_bam_files,
        out_file,
        chromosomes=default_chromosomes,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        min_normal_depth=0,
        min_tumour_depth=0,
        min_variant_depth=0,
        report_strand_counts=False,
        split_size=int(1e7),
        table_group=''):

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.TempOutputObj('regions_obj', 'regions'),
        value=biowrappers.components.variant_calling.utils.get_bam_regions(normal_bam_file, split_size, chromosomes=chromosomes)
    )

    tumour_input_files = {}

    for sample in tumour_bam_files:
        tumour_input_files[sample] = mgd.InputFile(tumour_bam_files[sample])

    workflow.transform(
        name='get_counts',
        axes=('regions',),
        ctx=med_ctx,
        func='biowrappers.components.snv_allele_counts.tasks.get_variant_position_counts',
        args=(
            mgd.InputFile(normal_bam_file),
            tumour_input_files,
            mgd.TempOutputFile('counts.h5', 'regions'),
            mgd.TempInputObj('regions_obj', 'regions'),
            table_group
        ),
        kwargs={
            'count_duplicates': count_duplicates,
            'min_bqual': min_bqual,
            'min_mqual': min_mqual,
            'min_normal_depth': min_normal_depth,
            'min_tumour_depth': min_tumour_depth,
            'min_variant_depth': min_variant_depth,
            'report_strand_counts': report_strand_counts
        }
    )

    workflow.transform(
        name='concatenate_counts',
        axes=(),
        ctx=med_ctx,
        func='biowrappers.components.io.hdf5.tasks.concatenate_tables',
        args=(
            mgd.TempInputFile('counts.h5', 'regions'),
            mgd.OutputFile(out_file),
        )
    )

    return workflow
