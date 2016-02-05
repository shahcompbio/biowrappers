from pypeliner.workflow import Workflow

import pypeliner
import biowrappers.components.io.bam.tasks as bam_tasks
import biowrappers.components.alignment.bwa.tasks as bwa_tasks

in_file = pypeliner.managed.InputFile('/home/andrew/Desktop/old/museq_test/input/DAH365B.chrom_21.bam')

ref_genome_fasta_file = pypeliner.managed.InputFile('/home/andrew/Desktop/old/museq_test/input/GRCh37-lite.fa')

tmp_dir = '/home/andrew/Desktop/realignment/tmp'

out_file = '/home/andrew/Desktop/realignment/aligned.bam'

read_1 = pypeliner.managed.TempFile('read_1', 'split')

read_2 = pypeliner.managed.TempFile('read_2', 'split')

read_1_sai = pypeliner.managed.TempFile('read_1.sai', 'split')

read_2_sai = pypeliner.managed.TempFile('read_2.sai', 'split')

workflow = Workflow()

workflow.transform(
    name='bam_to_fasta',
    axes=(),
    func=bam_tasks.split_bam_to_fastq, 
    args=(
        in_file,
        {
            1 : read_1.as_output(),
            2 : read_2.as_output(),
        },
        pypeliner.managed.TempSpace('bam_to_fastq'),
    ),
    kwargs={
        'split_size' : int(1e4)
    },
)

workflow.transform(
    name='aln_read_1',
    axes=('split',),
    func=bwa_tasks.run_aln,
    args=(
        read_1.as_input(),
        ref_genome_fasta_file,
        read_1_sai.as_output(),
    ),
)

workflow.transform(
    name='aln_read_2',
    axes=('split',),
    func=bwa_tasks.run_aln,
    args=(
        read_2.as_input(),
        ref_genome_fasta_file,
        read_2_sai.as_output(),
    ),
)

workflow.transform(
    name='sampe', 
    axes=('split',), 
    func=bwa_tasks.run_sampe, 
    args=(
        read_1.as_input(),
        read_2.as_input(),
        read_1_sai.as_input(),
        read_2_sai.as_input(),
        ref_genome_fasta_file,
        pypeliner.managed.TempOutputFile('aligned.bam', 'split'),
    ), 
)

workflow.transform(
    name='sort',
    axes=('split',),
    func=bam_tasks.sort_bam,
    args=(
        pypeliner.managed.TempInputFile('aligned.bam', 'split'),
        pypeliner.managed.TempOutputFile('sorted.bam', 'split')
    ),
)

workflow.transform(
    name='merge', 
    axes=(), 
    func=bam_tasks.merge_bams, 
    args=(
        pypeliner.managed.TempInputFile('sorted.bam', 'split'),
        pypeliner.managed.TempOutputFile('merged.bam')
    ),
)

workflow.transform(
    name='mark_duplicates',
    func=bam_tasks.mark_duplicates,
    args=(
        pypeliner.managed.TempInputFile('merged.bam'),
        pypeliner.managed.TempOutputFile('markdup.bam'),
    ),
)

workflow.commandline(
    name='copy',
    args=(
        'cp',
        pypeliner.managed.TempInputFile('markdup.bam'),
        pypeliner.managed.OutputFile(out_file),
    ), 
)

pyp = pypeliner.app.Pypeline([], {'maxjobs' : 6, 'submit' : 'local', 'tmpdir' : tmp_dir, 'nocleanup' : True, 'repopulate' : True})

pyp.run(workflow)