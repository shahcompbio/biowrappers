from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks

def vcf_mappability_annotation_pipeline(
    mappability_file,
    vcf_file,
    out_file,
    split_size=int(1e4),
    table_name='mappability'):

    workflow = Workflow()
    
    workflow.transform(
        name='split_vcf',
        ctx={'mem' : 1},
        func=vcf_tasks.split_vcf,
        args=(
            pypeliner.managed.InputFile(vcf_file),
            pypeliner.managed.TempOutputFile('split.vcf', 'split')
        ),
        kwargs={'lines_per_file' : split_size}
    )
    
    workflow.transform(
        name='annotate_db_status',
        axes=('split',),
        ctx={'mem' : 2},
        func=tasks.get_mappability,
        args=(
            pypeliner.managed.InputFile(mappability_file),
            pypeliner.managed.TempInputFile('split.vcf', 'split'),
            pypeliner.managed.TempOutputFile('mappability.h5', 'split'),
            table_name
        ),
    )
    
    workflow.transform(
        name='merge_tables',
        ctx={'mem' : 2},
        func=hdf5_tasks.concatenate_tables,
        args=(
            pypeliner.managed.TempInputFile('mappability.h5', 'split'),
            pypeliner.managed.OutputFile(out_file)
        )
    )
    
    return workflow
