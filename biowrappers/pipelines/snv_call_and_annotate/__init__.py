from pypeliner.workflow import Workflow

import os
import pypeliner

from biowrappers.components.utils import make_parent_directory
from biowrappers.components.variant_calling.utils import default_chromosomes

default_ctx = {'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2}

big_mem_ctx = {'mem': 8, 'num_retry': 3, 'mem_retry_increment': 2}


def call_and_annotate_pipeline(
        config,
        normal_bam_path,
        tumour_bam_paths,
        raw_data_dir,
        results_file,
        chromosomes=default_chromosomes):

    workflow = Workflow()

    workflow.setobj(pypeliner.managed.OutputChunks('tumour_sample_id', axes_origin=[0, ]), tumour_bam_paths.keys())

    variant_files = get_variant_files(chromosomes, config, raw_data_dir)

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
            func='biowrappers.components.variant_calling.nuseq.create_nuseq_classify_workflow',
            args=(
                normal_bam_file.as_input(),
                [pypeliner.managed.InputFile(x) for x in tumour_bam_paths.values()],
                ref_genome_fasta_file.as_input(),
                variant_files['snv']['vcf']['nuseq_multi_sample'].as_output()
            ),
            kwargs=config['nuseq_multi_sample']['kwargs']
        )

        workflow.transform(
            name='convert_nuseq_multi_sample_vcf_to_hdf5',
            axes=(),
            ctx=default_ctx,
            func="biowrappers.components.io.vcf.tasks.convert_vcf_to_hdf5",
            args=(
                variant_files['snv']['vcf']['nuseq_multi_sample'].as_input(),
                variant_files['snv']['hdf']['nuseq_multi_sample'].as_output(),
                '/snv/vcf/nuseq_multi_sample/all',
            ),
            kwargs={
                'score_callback': vcf_score_callbacks['snv']['nuseq']
            }
        )

    #===================================================================================================================
    # Single sample calling
    #===================================================================================================================
    if 'nuseq' in config:
        workflow.subworkflow(
            name='nuseq',
            axes=('tumour_sample_id',),
            func='biowrappers.components.variant_calling.nuseq.create_nuseq_classify_workflow',
            args=(
                normal_bam_file.as_input(),
                [tumour_bam_files.as_input(), ],
                ref_genome_fasta_file.as_input(),
                variant_files['snv']['vcf']['nuseq'].as_output()
            ),
            kwargs=config['nuseq']['kwargs']
        )

    if 'mutect' in config:
        workflow.subworkflow(
            name='mutect',
            axes=('tumour_sample_id',),
            func='biowrappers.components.variant_calling.mutect.create_mutect_workflow',
            args=(
                normal_bam_file.as_input(),
                tumour_bam_files.as_input(),
                ref_genome_fasta_file.as_input(),
                config['databases']['cosmic']['local_path'],
                config['databases']['dbsnp']['local_path'],
                variant_files['snv']['vcf']['mutect'].as_output()
            ),
            kwargs=config['mutect']['kwargs']
        )

    if 'strelka' in config:
        workflow.subworkflow(
            name='strelka',
            axes=('tumour_sample_id',),
            func='biowrappers.components.variant_calling.strelka.create_strelka_workflow',
            args=(
                normal_bam_file.as_input(),
                tumour_bam_files.as_input(),
                ref_genome_fasta_file.as_input(),
                variant_files['indel']['vcf']['strelka'].as_output(),
                variant_files['snv']['vcf']['strelka'].as_output()
            ),
            kwargs=config['strelka']['kwargs']
        )

    #===================================================================================================================
    # Convert vcf to hdf5
    #===================================================================================================================
    for var_type in variant_files:
        for prog in variant_files[var_type]['vcf']:
            if prog == 'nuseq_multi_sample':
                continue

            workflow.transform(
                name='convert_{0}_indel_{1}_to_hdf5'.format(prog, var_type),
                axes=('tumour_sample_id',),
                ctx=default_ctx,
                func="biowrappers.components.io.vcf.tasks.convert_vcf_to_hdf5",
                args=(
                    variant_files[var_type]['vcf'][prog].as_input(),
                    variant_files[var_type]['hdf'][prog].as_output(),
                    pypeliner.managed.Template(
                        '/{var_type}/vcf/{prog}/{{tumour_sample_id}}'.format(prog=prog, var_type=var_type),
                        'tumour_sample_id'
                    )
                ),
                kwargs={
                    'score_callback': vcf_score_callbacks[var_type][prog]
                }
            )

    #===================================================================================================================
    # Indel annotation
    #===================================================================================================================
    workflow.transform(
        name='merge_indels',
        ctx=big_mem_ctx,
        func='biowrappers.components.io.vcf.tasks.vcf_tasks.merge_vcfs',
        args=(
            [x.as_input() for x in variant_files['indel']['vcf'].values()],
            pypeliner.managed.TempOutputFile('all.indel.vcf')
        )
    )

    workflow.transform(
        name='finalise_indels',
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
        args=(
            pypeliner.managed.TempInputFile('all.indel.vcf'),
            pypeliner.managed.TempOutputFile('all.indel.vcf.gz')
        )
    )

    workflow.subworkflow(
        name='annotate_indels',
        axes=(),
        func=create_annotation_workflow,
        args=(
            config,
            pypeliner.managed.TempInputFile('all.indel.vcf.gz'),
            pypeliner.managed.TempOutputFile('indel_annotations.h5'),
            os.path.join(raw_data_dir, 'indel'),
        ),
        kwargs={
            'variant_type': 'indel'
        }
    )

    #===================================================================================================================
    # SNV
    #===================================================================================================================
    workflow.transform(
        name='merge_snvs',
        ctx=big_mem_ctx,
        func="biowrappers.components.io.vcf.tasks.merge_vcfs",
        args=(
            [x.as_input() for x in variant_files['snv']['vcf'].values()],
            pypeliner.managed.TempOutputFile('all.snv.vcf')
        )
    )

    workflow.transform(
        name='finalise_snvs',
        func="biowrappers.components.io.vcf.tasks.finalise_vcf",
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
            os.path.join(raw_data_dir, 'snv'),
        ),
        kwargs={
            'variant_type': 'snv'
        }
    )

    workflow.subworkflow(
        name='normal_snv_counts',
        func='biowrappers.components.variant_calling.snv_allele_counts.create_snv_allele_counts_for_vcf_targets_workflow',
        args=(
            normal_bam_file.as_input(),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.OutputFile(os.path.join(raw_data_dir, 'snv', 'counts', 'normal.h5')),
        ),
        kwargs=get_kwargs(config['snv_counts']['kwargs'], '/snv/counts/normal')
    )

    workflow.subworkflow(
        name='tumour_snv_counts',
        axes=('tumour_sample_id',),
        func='biowrappers.components.variant_calling.snv_allele_counts.create_snv_allele_counts_for_vcf_targets_workflow',
        args=(
            tumour_bam_files.as_input(),
            pypeliner.managed.TempInputFile('all.snv.vcf.gz'),
            pypeliner.managed.OutputFile(
                os.path.join(raw_data_dir, 'snv', 'counts', '{tumour_sample_id}.h5'), 'tumour_sample_id')
        ),
        kwargs=get_kwargs(
            config['snv_counts']['kwargs'],
            pypeliner.managed.Template('/snv/counts/{tumour_sample_id}', 'tumour_sample_id')
        )
    )

    #===================================================================================================================
    # Create final output
    #===================================================================================================================
    tables = [
        pypeliner.managed.TempInputFile('indel_annotations.h5'),
        pypeliner.managed.TempInputFile('snv_annotations.h5'),
        pypeliner.managed.InputFile(os.path.join(raw_data_dir, 'snv', 'counts', 'normal.h5')),
        pypeliner.managed.InputFile(
            os.path.join(raw_data_dir, 'snv', 'counts', '{tumour_sample_id}.h5'), 'tumour_sample_id'),
    ]

    for var_type in variant_files:
        for prog in variant_files[var_type]['hdf']:
            tables.append(variant_files[var_type]['hdf'][prog].as_input())

    workflow.transform(
        name='build_results_file',
        ctx=default_ctx,
        func='biowrappers.components.io.hdf5.tasks.concatenate_tables',
        args=(
            tables,
            pypeliner.managed.OutputFile(results_file)
        ),
        kwargs={
            'drop_duplicates': True,
        }
    )

    return workflow


def create_annotation_workflow(
        config,
        in_vcf_file,
        out_file,
        raw_data_dir,
        variant_type='snv',
        docker_config={},
        snpeff_docker={},
        vcftools_docker={},
):

    annotators = (
        'cosmic_status',
        'dbsnp_status',
        'mappability',
        'snpeff',
        'tri_nucleotide_context'
    )

    result_files = {}

    kwargs = {}

    for a in annotators:
        kwargs[a] = get_kwargs(config[a]['kwargs'], '/{0}/{1}'.format(variant_type, a))

        result_files[a] = pypeliner.managed.File(os.path.join(raw_data_dir, '{0}.csv.gz'.format(a)))

    if not os.path.isdir(raw_data_dir):
        os.mkdir(raw_data_dir)

    assert os.path.isdir(raw_data_dir)

    workflow = Workflow()

    workflow.subworkflow(
        name='cosmic_status',
        func='biowrappers.components.variant_calling.annotated_db_status.create_vcf_db_annotation_workflow',
        ctx=dict(mem=4, mem_retry_increment=2, **docker_config),
        args=(
            config['databases']['cosmic']['local_path'],
            pypeliner.managed.InputFile(in_vcf_file),
            result_files['cosmic_status'].as_output(),
        ),
        kwargs=config["cosmic_status"]['kwargs']
    )

    workflow.subworkflow(
        name='dbsnp_status',
        func='biowrappers.components.variant_calling.annotated_db_status.create_vcf_db_annotation_workflow',
        ctx=dict(mem=4, mem_retry_increment=2, **docker_config),
        args=(
            config['databases']['dbsnp']['local_path'],
            pypeliner.managed.InputFile(in_vcf_file),
            result_files['dbsnp_status'].as_output(),
        ),
        kwargs=config["dbsnp_status"]['kwargs']
    )

    workflow.subworkflow(
        name='mappability',
        func='biowrappers.components.variant_calling.mappability.create_vcf_mappability_annotation_workflow',
        ctx=dict(mem=4, mem_retry_increment=2, **docker_config),
        args=(
            config['databases']['mappability']['local_path'],
            pypeliner.managed.InputFile(in_vcf_file, extensions=['.tbi']),
            result_files['mappability'].as_output(),
        ),
        kwargs=config["mappability"]['kwargs']
    )

    workflow.subworkflow(
        name='snpeff',
        func='biowrappers.components.variant_calling.snpeff.create_snpeff_annotation_workflow',
        ctx=dict(mem=4, mem_retry_increment=2, **docker_config),
        args=(
            config['databases']['snpeff']['db'],
            config['databases']['snpeff']['data_dir'],
            pypeliner.managed.InputFile(in_vcf_file),
            result_files['snpeff'].as_output(),
        ),

        kwargs=dict(snpeff_docker=snpeff_docker,  **kwargs['snpeff'])
    )

    workflow.subworkflow(
        name='tri_nucleotide_context',
        func='biowrappers.components.variant_calling.tri_nucleotide_context.create_vcf_tric_nucleotide_annotation_workflow',
        ctx=dict(mem=4, mem_retry_increment=2, **docker_config),
        args=(
            config['databases']['ref_genome']['local_path'],
            pypeliner.managed.InputFile(in_vcf_file),
            result_files['tri_nucleotide_context'].as_output(),
        ),
        kwargs=config["tri_nucleotide_context"]['kwargs']
    )

    workflow.transform(
        name='build_results_file',
        ctx=dict(mem=4, mem_retry_increment=2, **docker_config),
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            [x.as_input() for x in result_files.values()],
            pypeliner.managed.OutputFile(out_file, extensions=[".yaml"]),
        )
    )

    return workflow


def get_variant_files(chromosomes, config, raw_data_dir):
    indel_hdf_files = {}

    indel_vcf_files = {}

    snv_hdf_files = {}

    snv_vcf_files = {}

    indel_progs = ('strelka', )

    indel_progs = [x for x in indel_progs if x in config]

    snv_progs = ('nuseq', 'mutect', 'strelka', )

    snv_progs = [x for x in snv_progs if x in config]

    for prog in snv_progs:
        config[prog]['kwargs']['chromosomes'] = chromosomes

        snv_hdf_files[prog] = pypeliner.managed.File(
            get_sample_out_file(prog, 'h5', raw_data_dir, 'snv'),
            'tumour_sample_id'
        )

        snv_vcf_files[prog] = pypeliner.managed.File(
            get_sample_out_file(prog, 'vcf.gz', raw_data_dir, 'snv'),
            'tumour_sample_id'
        )

        if prog in indel_progs:
            indel_hdf_files[prog] = pypeliner.managed.File(
                get_sample_out_file(prog, 'h5', raw_data_dir, 'indel'),
                'tumour_sample_id'
            )

            indel_vcf_files[prog] = pypeliner.managed.File(
                get_sample_out_file(prog, 'vcf.gz', raw_data_dir, 'indel'),
                'tumour_sample_id'
            )

    if 'nuseq_multi_sample' in config:
        config['nuseq_multi_sample']['kwargs']['chromosomes'] = chromosomes

        snv_hdf_files['nuseq_multi_sample'] = pypeliner.managed.File(
            get_sample_out_file('nuseq_multi_sample', 'h5', raw_data_dir, 'snv').format(tumour_sample_id='all')
        )

        snv_vcf_files['nuseq_multi_sample'] = pypeliner.managed.File(
            get_sample_out_file('nuseq_multi_sample', 'vcf.gz', raw_data_dir, 'snv').format(tumour_sample_id='all')
        )

    return {
        'indel': {
            'hdf': indel_hdf_files,
            'vcf': indel_vcf_files,
        },
        'snv': {
            'hdf': snv_hdf_files,
            'vcf': snv_vcf_files
        }
    }


def get_sample_out_file(cmd, ext, out_dir, variant_type):

    out_file = os.path.join(out_dir, variant_type, cmd, '{{tumour_sample_id}}.{0}'.format(ext))

    make_parent_directory(out_file)

    return out_file


def get_kwargs(config, table_name):

    config = config.copy()

    config['table_name'] = table_name

    return config


def nuseq_callback(record):
    return record.INFO['PS']


def strelka_indel_callback(record):
    return record.INFO['QSI']


def strelka_snv_callback(record):
    return record.INFO['QSS']

vcf_score_callbacks = {
    'indel': {
        'strelka': strelka_indel_callback,
    },
    'snv': {
        'mutect': None,
        'nuseq': nuseq_callback,
        'strelka': strelka_snv_callback,
    }
}
