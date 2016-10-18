from pypeliner.workflow import Workflow

import pypeliner

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.variant_calling.utils as utils
import tasks

def create_vcf_mappability_annotation_workflow(
    mappability_file,
    vcf_file,
    out_file,
    chromosomes=default_chromosomes,
    split_size=int(1e7),
    table_name='mappability'):

    workflow = Workflow()
    
    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('regions_obj', 'regions'),
        value=utils.get_vcf_regions(vcf_file, split_size, chromosomes=chromosomes, zero_based=True, half_open=True)
    )
    
    workflow.transform(
        name='annotate_db_status',
        axes=('split',),
        ctx={'mem' : 2, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=tasks.get_mappability,
        args=(
            pypeliner.managed.InputFile(mappability_file),
            pypeliner.managed.TempInputFile('split.vcf', 'split'),
            pypeliner.managed.TempOutputFile('mappability.h5', 'split'),
            table_name
        ),
        kwargs={
            'region': pypeliner.managed.TempInputObj('regions_obj', 'regions'),
        },
    )
    
    workflow.transform(
        name='merge_tables',
        ctx={'mem' : 2, 'num_retry' : 3, 'mem_retry_increment' : 2},
        func=hdf5_tasks.concatenate_tables,
        args=(
            pypeliner.managed.TempInputFile('mappability.h5', 'split'),
            pypeliner.managed.OutputFile(out_file)
        )
    )
    
    return workflow
