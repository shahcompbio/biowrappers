from pypeliner.workflow import Workflow

import pypeliner.managed as mgd

import tasks


def create_tophat_transcriptome_index_workflow(
    ref_genome_fasta_file,
    transcript_gtf_file,
    ref_genome_index_prefix,
    transcriptome_index_prefix,
    copy_ref_genome=False
):

    workflow = Workflow()
    local_ref_genome_fasta_path = ref_genome_index_prefix + '.fa'
    if copy_ref_genome:
        workflow.commandline(
            name='copy_genome',
            ctx={'local': True},
            args=(
                'cp',
                mgd.InputFile(ref_genome_fasta_file),
                mgd.OutputFile(local_ref_genome_fasta_path),
            ),
        )
    else:
        workflow.commandline(
            name='link_genome',
            ctx={'local': True},
            args=(
                'ln',
                '-s',
                mgd.InputFile(ref_genome_fasta_file),
                mgd.OutputFile(local_ref_genome_fasta_path),
            ),
        )
    workflow.transform(
        name='build_bowtie_index',
        ctx={'mem': 8, 'num_retry': 3, 'mem_retry_increment': 8},
        func=tasks.build_genome_index,
        args=(
            mgd.InputFile(local_ref_genome_fasta_path),
            mgd.OutputFile(ref_genome_index_prefix),
        )
    )
    workflow.transform(
        name='build_tophat_index',
        ctx={'mem': 8, 'num_retry': 3, 'mem_retry_increment': 8},
        func=tasks.build_transcriptome_index,
        args=(
            mgd.InputFile(ref_genome_index_prefix),
            mgd.InputFile(transcript_gtf_file),
            mgd.OutputFile(transcriptome_index_prefix),
        )
    )
    return workflow
