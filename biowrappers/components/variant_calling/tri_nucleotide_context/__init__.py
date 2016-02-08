from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.io.vcf.tasks as vcf_tasks
import tasks

def vcf_tric_nucleotide_annotation_pipeline(
    ref_genome_fasta_file,
    vcf_file,
    out_file,
    split_size=int(1e4),
    table_name='tri_nucleotide_context'):
    
    workflow = Workflow()
    
    workflow.transform(
        name='split_vcf',
        ctx={'mem' : 2},
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
        ctx={'mem' : 4},
        func=tasks.get_tri_nucelotide_context,
        args=(
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.TempInputFile('split.vcf', 'split'),
            pypeliner.managed.TempOutputFile('tri_nucleotide_context.h5', 'split'),
            table_name
        )
    )
    
    workflow.transform(
        name='merge_tables',
        ctx={'mem' : 2},
        func=hdf5_tasks.concatenate_tables,
        args=(
            pypeliner.managed.TempInputFile('tri_nucleotide_context.h5', 'split'),
            pypeliner.managed.OutputFile(out_file)
        )
    )
    
    return workflow
