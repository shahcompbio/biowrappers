import pypeliner
import pypeliner.managed as mgd


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def create_vcf_mappability_annotation_workflow(
        mappability_file,
        vcf_file,
        out_file,
        chromosomes=default_chromosomes,
        split_size=int(1e7),
):

    ctx = {'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2}

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='get_regions',
        ret=mgd.TempOutputObj('regions_obj', 'regions'),
        ctx=ctx,
        func='biowrappers.components.variant_calling.utils.get_vcf_regions',
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
        func='biowrappers.components.variant_calling.mappability.tasks.get_mappability',
        args=(
            mappability_file,
            mgd.InputFile(vcf_file, extensions=['.tbi']),
            mgd.TempOutputFile('mappability.csv.gz', 'regions')
        ),
        kwargs={
            'region': mgd.TempInputObj('regions_obj', 'regions'),
        },
    )

    workflow.transform(
        name='merge_tables',
        ctx=ctx,
        func='biowrappers.components.io.csv.tasks.concatenate_csv',
        args=(
            mgd.TempInputFile('mappability.csv.gz', 'regions'),
            mgd.OutputFile(out_file)
        )
    )

    return workflow
