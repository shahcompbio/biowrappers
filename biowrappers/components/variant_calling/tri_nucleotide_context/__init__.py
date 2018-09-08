import pypeliner
import pypeliner.managed as mgd

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks


def create_vcf_tric_nucleotide_annotation_workflow(
        ref_genome_fasta_file,
        vcf_file,
        out_file,
        docker_config=None,
        hdf5_output=True,
        split_size=int(1e4),
        table_name='tri_nucleotide_context'):

    ctx = {'num_retry': 3, 'mem_retry_increment': 2}
    if docker_config:
        ctx.update(docker_config)

    if hdf5_output:
        merged_file = mgd.File(out_file)
    else:
        merged_file = mgd.TempFile('merged.h5')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='split_vcf',
        ctx=dict(mem=2, **ctx),
        func=vcf_tasks.split_vcf,
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
        func=tasks.get_tri_nucelotide_context,
        args=(
            ref_genome_fasta_file,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('tri_nucleotide_context.h5', 'split'),
            table_name
        )
    )

    workflow.transform(
        name='merge_tables',
        ctx=dict(mem=2, **ctx),
        func=hdf5_tasks.concatenate_tables,
        args=(
            mgd.TempInputFile('tri_nucleotide_context.h5', 'split'),
            merged_file.as_output()
        )
    )

    if not hdf5_output:
        workflow.transform(
            name='convert_to_tsv',
            ctx=dict(mem=2, **ctx),
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
