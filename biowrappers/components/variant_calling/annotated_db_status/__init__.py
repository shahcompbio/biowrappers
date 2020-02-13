import pypeliner
import pypeliner.managed as mgd

def create_vcf_db_annotation_workflow(
        db_vcf_file,
        target_vcf_file,
        out_file,
        docker_config={},
        split_size=int(1e4)):

    ctx = dict(mem=2, num_retry=3, mem_retry_increment=2, **docker_config)

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='split_vcf',
        ctx=ctx,
        func='biowrappers.components.io.vcf.tasks.split_vcf',
        args=(
            mgd.InputFile(target_vcf_file),
            mgd.TempOutputFile('split.vcf', 'split')
        ),
        kwargs={'lines_per_file': split_size}
    )

    workflow.transform(
        name='annotate_db_status',
        axes=('split',),
        ctx=ctx,
        func='biowrappers.components.variant_calling.annotated_db_status.tasks.annotate_db_status',
        args=(
            db_vcf_file,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('annotated.csv.gz', 'split',
                               extensions=['.yaml'])
        )
    )

    workflow.transform(
        name='merge_tables',
        ctx=ctx,
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('annotated.csv.gz', 'split'),
            mgd.OutputFile(out_file, extensions=['.yaml'])
        )
    )

<<<<<<< HEAD
    if not hdf5_output:
        workflow.transform(
            name='convert_to_tsv',
            ctx=ctx,
            func='single_cell.utils.hdfutils.convert_hdf_to_csv',
            args=(
                merged_file.as_input(),
                mgd.OutputFile(out_file),
            ),
            kwargs={
                'compress': True,
            }
        )

=======
>>>>>>> 3345feb... stripped out hd5 and adde csv for vcf annotation
    return workflow
