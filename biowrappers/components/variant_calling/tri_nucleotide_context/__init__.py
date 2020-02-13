import pypeliner
import pypeliner.managed as mgd


def create_vcf_tric_nucleotide_annotation_workflow(
        ref_genome_fasta_file,
        vcf_file,
        out_file,
        docker_config=None,
        split_size=int(1e4),
        table_name='tri_nucleotide_context'):

    ctx = {'num_retry': 3, 'mem_retry_increment': 2}
    if docker_config:
        ctx.update(docker_config)


    merged_file = mgd.TempFile('merged.csv.gz')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='split_vcf',
        ctx=dict(mem=2, **ctx),
        func='biowrappers.components.io.vcf.tasks.split_vcf',
        args=(
            mgd.InputFile(vcf_file),
            mgd.TempOutputFile('split.vcf', 'split')
        ),
        kwargs={'lines_per_file': split_size}
    )

    workflow.transform(
        name='annotate_db_status',
        axes=('split',),
        ctx=dict(mem=4, **ctx),
        func='biowrappers.components.variant_calling.tri_nucleotide_context.tasks.get_tri_nucelotide_context',
        args=(
            ref_genome_fasta_file,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('tri_nucleotide_context.csv.gz', 'split',
                               extensions=['.yaml']),
            table_name
        )
    )

    workflow.transform(
        name='merge_tables',
        ctx=dict(mem=2, **ctx),
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('tri_nucleotide_context.csv.gz', 'split'),
            mgd.OutputFile(out_file, extensions=['.yaml']))
    )


    return workflow
