from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks

def snpeff_pipeline(
    target_vcf_file,
    out_file,
    data_base='GRCh37.75',
    memory=4,
    split_size=int(1e3),
    table_name='snpeff'):
    
    workflow = Workflow()
    
    workflow.transform(
        name='split_vcf',
        ctx={'mem' : 1},
        func=vcf_tasks.split_vcf,
        args=(
            pypeliner.managed.InputFile(target_vcf_file),
            pypeliner.managed.TempOutputFile('split.vcf', 'split')
        ),
        kwargs={'lines_per_file' : split_size}
    )
    
    workflow.transform(
        name='run_snpeff',
        axes=('split',),
        ctx={'mem' : memory + 2},
        func=tasks.run_snpeff,
        args=(
            pypeliner.managed.TempInputFile('split.vcf', 'split'),
            pypeliner.managed.TempOutputFile('snpeff.vcf', 'split')
        ),
        kwargs={'data_base' : data_base, 'memory' : memory}
    )
    
    workflow.transform(
        name='convert_vcf_to_table',
        axes=('split',),
        ctx={'mem' : 2},
        func=tasks.convert_vcf_to_table,
        args=(
            pypeliner.managed.TempInputFile('snpeff.vcf', 'split'),
            pypeliner.managed.TempOutputFile('snpeff.h5', 'split'),
            table_name
        )
    )
    
    workflow.transform(
        name='concatenate_tables',
        ctx={'mem' : 2},
        func=hdf5_tasks.concatenate_tables,
        args=(
            pypeliner.managed.TempInputFile('snpeff.h5', 'split'),
            pypeliner.managed.OutputFile(out_file)
        )
    )
    
    return workflow