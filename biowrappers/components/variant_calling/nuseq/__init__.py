from pypeliner.workflow import Workflow

import pypeliner

import biowrappers.components.io.vcf.tasks as vcf_tasks
import biowrappers.components.variant_calling.utils as utils

import tasks


def create_nuseq_classify_workflow(
        normal_bam_file,
        tumour_bam_files,
        ref_genome_fasta_file,
        out_file,
        chromosomes=utils.default_chromosomes,
        chunk_size=int(1e5),
        indel_threshold=0.05,
        min_normal_depth=1,
        min_tumour_depth=1,
        min_somatic_probability=0.0,
        split_size=int(1e7)):

    workflow = Workflow()

    workflow.transform(
        name='get_regions',
        func=utils.get_bam_regions,
        ret=pypeliner.managed.TempOutputObj('regions', 'regions'),
        args=(
            pypeliner.managed.InputFile(tumour_bam_files[0]),
            split_size,
        ),
        kwargs={
            'chromosomes': chromosomes,
        },
    )

    workflow.transform(
        name='run_classify',
        axes=('regions',),
        ctx={'mem': 6, 'num_retry': 3, 'mem_retry_increment': 2, 'io': 1},
        func=tasks.run_classify,
        args=(
            pypeliner.managed.InputFile(normal_bam_file),
            [pypeliner.managed.InputFile(x) for x in tumour_bam_files],
            pypeliner.managed.InputFile(ref_genome_fasta_file),
            pypeliner.managed.TempInputObj('regions', 'regions'),
            pypeliner.managed.TempOutputFile('classified.h5', 'regions')
        ),
        kwargs={
            'chunk_size': chunk_size,
            'min_normal_depth': min_normal_depth,
            'min_tumour_depth': min_tumour_depth,
            'min_somatic_probability': min_somatic_probability
        }
    )

    workflow.transform(
        name='write_vcf',
        axes=('regions',),
        ctx={'mem': 4, 'num_retry': 3, 'mem_retry_increment': 2},
        func=tasks.write_vcf,
        args=(
            pypeliner.managed.TempInputFile('classified.h5', 'regions'),
            pypeliner.managed.TempOutputFile('classified.vcf', 'regions')
        ),
        kwargs={'indel_threshold': indel_threshold}
    )

    workflow.transform(
        name='merge_vcf',
        ctx={'mem': 8, 'num_retry': 3, 'mem_retry_increment': 2},
        func=vcf_tasks.concatenate_vcf,
        args=(
            pypeliner.managed.TempInputFile('classified.vcf', 'regions'),
            pypeliner.managed.TempOutputFile('merged.vcf')
        )
    )

    workflow.transform(
        name='filter_snvs',
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2},
        func=vcf_tasks.filter_vcf,
        args=(
            pypeliner.managed.TempInputFile('merged.vcf'),
            pypeliner.managed.TempOutputFile('merged.filtered.vcf')
        )
    )

    workflow.transform(
        name='finalise',
        func=vcf_tasks.finalise_vcf,
        args=(
            pypeliner.managed.TempInputFile('merged.filtered.vcf'),
            pypeliner.managed.OutputFile(out_file)
        )
    )

    return workflow
