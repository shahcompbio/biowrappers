from pypeliner.workflow import Workflow

import pypeliner.managed as mgd


def create_snpeff_annotation_workflow(
        db,
        data_dir,
        target_vcf_file,
        out_file,
        base_docker={},
        snpeff_docker={},
        classic_mode=True,
        split_size=int(1e3),
        table_name='snpeff'):

    ctx = {'num_retry': 3, 'mem_retry_increment': 2}

    if base_docker:
        ctx.update(base_docker)

    workflow = Workflow()

    workflow.transform(
        name='split_vcf',
        ctx=dict(mem=2, **ctx),
        func='biowrappers.components.io.vcf.tasks.split_vcf',
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
        func='biowrappers.components.variant_calling.snpeff.tasks.run_snpeff',
        args=(
            db,
            data_dir,
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('snpeff.vcf', 'split')
        ),
        kwargs={
            'classic_mode': classic_mode,
            'docker_config': snpeff_docker
        }
    )

    workflow.transform(
        name='convert_vcf_to_csv',
        axes=('split',),
        ctx=dict(mem=4, **ctx),
        func='biowrappers.components.variant_calling.snpeff.tasks.convert_vcf_to_table',
        args=(
            mgd.TempInputFile('snpeff.vcf', 'split'),
            mgd.TempOutputFile('snpeff.csv.gz', 'split', extensions=['.yaml']),
            table_name
        )
    )

    workflow.transform(
        name='concatenate_tables',
        ctx=dict(mem=4, **ctx),
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('snpeff.csv.gz', 'split'),
            mgd.OutputFile(out_file, extensions=['.yaml'])
        )
    )

    return workflow