import pypeliner
import pypeliner.managed as mgd

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks


def create_vcf_db_annotation_workflow(
        db_vcf_file,
        target_vcf_file,
        out_file,
        table_name,
        hdf5_output=True,
        split_size=int(1e4)):

    if hdf5_output:
        merged_file = mgd.File(out_file)

    else:
        merged_file = mgd.TempFile('merged.h5')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='split_vcf',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=vcf_tasks.split_vcf,
        args=(
            mgd.InputFile(target_vcf_file),
            mgd.TempOutputFile('split.vcf', 'split')
        ),
        kwargs={'lines_per_file': split_size}
    )

    workflow.transform(
        name='annotate_db_status',
        axes=('split',),
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=tasks.annotate_db_status,
        args=(
            db_vcf_file,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('annotated.h5', 'split'),
            table_name
        )
    )

    workflow.transform(
        name='merge_tables',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=hdf5_tasks.concatenate_tables,
        args=(
            mgd.TempInputFile('annotated.h5', 'split'),
            merged_file.as_output()
        )
    )

    if not hdf5_output:
        workflow.transform(
            name='convert_to_tsv',
            ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
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
