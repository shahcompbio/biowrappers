import pypeliner
import pypeliner.managed as mgd

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.variant_calling.utils as utils
import tasks


def create_vcf_mappability_annotation_workflow(
        mappability_file,
        vcf_file,
        out_file,
        docker_config={},
        chromosomes=default_chromosomes,
        hdf5_output=True,
        split_size=int(1e7),
        table_name='mappability',
        base_docker={}):

    ctx = {'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2}
    if docker_config:
        ctx.update(docker_config)

    if hdf5_output:
        merged_file = mgd.File(out_file)
    else:
        merged_file = mgd.TempFile('merged.h5')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='get_regions',
        ret=mgd.TempOutputObj('regions_obj', 'regions'),
        ctx=ctx,
        func=utils.get_vcf_regions,
        args=(
            mgd.InputFile(vcf_file, extensions=['.tbi']),
            split_size,
        ),
        kwargs={
            'chromosomes': chromosomes,
        },
    )

    workflow.transform(
        name='annotate_db_status',
        axes=('regions',),
        ctx=ctx,
        func=tasks.get_mappability,
        args=(
            mappability_file,
            mgd.InputFile(vcf_file, extensions=['.tbi']),
            mgd.TempOutputFile('mappability.h5', 'regions'),
            table_name
        ),
        kwargs={
            'region': mgd.TempInputObj('regions_obj', 'regions'),
        },
    )

    workflow.transform(
        name='merge_tables',
        ctx=ctx,
        func=hdf5_tasks.concatenate_tables,
        args=(
            mgd.TempInputFile('mappability.h5', 'regions'),
            merged_file.as_output()
        )
    )

    if not hdf5_output:
        workflow.transform(
            name='convert_to_tsv',
            ctx=ctx,
            func=hdf5_tasks.convert_hdf5_to_tsv,
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
