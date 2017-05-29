from pypeliner.workflow import Workflow

import pypeliner
import pypeliner.managed as mgd

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks


def create_snpeff_annotation_workflow(
        db,
        target_vcf_file,
        out_file,
        hdf5_output=True,
        split_size=int(1e3),
        table_name='snpeff'):

    workflow = Workflow()

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
        name='run_snpeff',
        axes=('split',),
        ctx={'mem': 8, 'num_retry': 3, 'mem_retry_increment': 2},
        func=tasks.run_snpeff,
        args=(
            db,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('snpeff.vcf', 'split')
        )
    )

    if hdf5_output:
        workflow.transform(
            name='convert_vcf_to_table',
            axes=('split',),
            ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
            func=tasks.convert_vcf_to_table,
            args=(
                mgd.TempInputFile('snpeff.vcf', 'split'),
                mgd.TempOutputFile('snpeff.h5', 'split'),
                table_name
            )
        )

        workflow.transform(
            name='concatenate_tables',
            ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
            func=hdf5_tasks.concatenate_tables,
            args=(
                mgd.TempInputFile('snpeff.h5', 'split'),
                mgd.OutputFile(out_file)
            )
        )

    else:
        workflow.transform(
            name='run_snpeff',
            axes=('split',),
            ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
            func=vcf_tasks.finalise_vcf,
            args=(
                mgd.TempInputFile('snpeff.vcf', 'split'),
                mgd.TempOutputFile('snpeff.vcf.gz', 'split'),
            )
        )

        workflow.transform(
            name='merge_vcf',
            axes=(),
            ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
            func=vcf_tasks.concatenate_vcf,
            args=(
                mgd.TempInputFile('snpeff.vcf.gz', 'split'),
                mgd.OutputFile(out_file),
            )
        )

    return workflow
