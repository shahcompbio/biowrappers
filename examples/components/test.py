from pypeliner.workflow import Workflow

import pypeliner
import biowrappers.components.io.bam.tasks as tasks

file_name = '/home/andrew/Desktop/old/museq_test/input/DAH365B.chrom_21.bam'

read_1 = '/home/andrew/Desktop/realignment/read_1.{split}.fastq.gz'

read_2 = '/home/andrew/Desktop/realignment/read_2.{split}.fastq.gz'

workflow = Workflow()

workflow.transform(
    name='bam_to_fasta',
    axes=(),
    func=tasks.split_bam_to_fastq, 
    args=(
        pypeliner.managed.InputFile(file_name),
        {
            1 : pypeliner.managed.OutputFile(read_1, 'split'),
            2 : pypeliner.managed.OutputFile(read_2, 'split')
        },
        pypeliner.managed.TempSpace('bam_to_fastq'),
    ),
    kwargs={
        'split_size' : int(1e4)
    }
)

pyp = pypeliner.app.Pypeline([], {'submit' : 'local'})

pyp.run(workflow)