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
        base_docker={},
        snpeff_docker={},
        vcftools_docker={},
        classic_mode=True,
        hdf5_output=True,
        split_size=int(1e3),
        table_name='snpeff'):

    ctx = {'num_retry': 3, 'mem_retry_increment': 2}

    if base_docker:
        ctx.update(base_docker)

    workflow = Workflow()

    workflow.transform(
        name='split_vcf',
        ctx=dict(mem=2, **ctx),
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
        ctx=dict(mem=8, **ctx),
        func=tasks.run_snpeff,
        args=(
            db,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('snpeff.vcf', 'split')
        ),
        kwargs={
            'classic_mode': classic_mode,
            'docker_config': snpeff_docker
        }
    )

    if hdf5_output:
        workflow.transform(
            name='convert_vcf_to_table',
            axes=('split',),
            ctx=dict(mem=4, **ctx),
            func=tasks.convert_vcf_to_table,
            args=(
                mgd.TempInputFile('snpeff.vcf', 'split'),
                mgd.TempOutputFile('snpeff.h5', 'split'),
                table_name
            )
        )

        workflow.transform(
            name='concatenate_tables',
            ctx=dict(mem=4, **ctx),
            func=hdf5_tasks.concatenate_tables,
            args=(
                mgd.TempInputFile('snpeff.h5', 'split'),
                mgd.OutputFile(out_file)
            )
        )

    else:
        workflow.transform(
            name='compress_split_vcf',
            axes=('split',),
            ctx=dict(mem=2, **ctx),
            func=vcf_tasks.finalise_vcf,
            args=(
                mgd.TempInputFile('snpeff.vcf', 'split'),
                mgd.TempOutputFile('snpeff.vcf.gz', 'split', extensions=['.tbi', '.csi']),
            ),
            kwargs={'docker_config':vcftools_docker}
        )

        workflow.transform(
            name='merge_vcf',
            axes=(),
            ctx=dict(mem=4, **ctx),
            func=vcf_tasks.concatenate_vcf,
            args=(
                mgd.TempInputFile('snpeff.vcf.gz', 'split', extensions=['.tbi', '.csi']),
                mgd.OutputFile(out_file),
            ),
            kwargs={'docker_config': vcftools_docker,
                    'allow_overlap': True}
        )

    return workflow
