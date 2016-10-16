from biowrappers.components.io.download import create_download_workflow
from pypeliner.workflow import Workflow

import os
import pypeliner.managed as mgd

import biowrappers.components.rna.kallisto.tasks
import biowrappers.components.rna.salmon.tasks
import biowrappers.components.rna.star.tasks
import biowrappers.components.rna.tophat.workflow


def download_external_files(config):
    download_keys = [x for x in config if 'url' in config[x]]
    urls = dict(zip(
        download_keys,
        [config[x]['url'] for x in download_keys],
    ))
    print urls
    downloaded_files = dict(zip(
        urls.keys(),
        [config[x]['local_path'] for x in urls.keys()],
    ))

    workflow = Workflow()
    workflow.setobj(
        obj=mgd.TempOutputObj('url', 'files'),
        value=urls,
    )
    workflow.subworkflow(
        name='download',
        func=create_download_workflow,
        axes=('files',),
        args=(
            mgd.TempInputObj('url', 'files'),
            mgd.TempOutputFile('download.gz', 'files'),
        ),
    )
    workflow.commandline(
        name='unzip',
        axes=('files',),
        args=(
            'gzip',
            '-cd',
            mgd.TempInputFile('download.gz', 'files'),
            '>',
            mgd.OutputFile('unzipped', 'files', fnames=downloaded_files)
        ),
    )
    return workflow


def build_indexes(config):
    workflow = Workflow()
    if 'kallisto' in config:
        workflow.transform(
            name='build_kallisto_index',
            func=biowrappers.components.rna.kallisto.tasks.build_index,
            ctx={'mem': 32, 'num_retry': 3, 'mem_retry_increment': 8},
            args=(
                mgd.OutputFile(config['kallisto']['index']),
                mgd.InputFile(config['transcriptome_fasta_file']['local_path']),
            ),
            kwargs={
                'kmer_length': config['kallisto']['kmer_length']
            }
        )
    if 'salmon' in config:
        workflow.transform(
            name='build_salmon_index',
            func=biowrappers.components.rna.salmon.tasks.build_index,
            ctx={'mem': 32, 'num_retry': 3, 'mem_retry_increment': 8},
            args=(
                mgd.OutputFile(os.path.join(config['salmon']['index'], 'index.finished')),
                mgd.InputFile(config['transcriptome_fasta_file']['local_path']),
            ),
            kwargs={
                'kmer_length': config['salmon']['kmer_length'],
                'gencode': config['salmon'].get('gencode', False),
            }
        )
    if 'star' in config:
        workflow.transform(
            name='build_star_index',
            func=biowrappers.components.rna.star.tasks.build_index,
            ctx={'mem': 32, 'num_retry': 3, 'mem_retry_increment': 8, 'local': config['star'].get('local', False)},
            args=(
                mgd.OutputFile(os.path.join(config['star']['index'], 'index.finished')),
                mgd.InputFile(config['ref_genome_fasta_file']['local_path']),
                mgd.InputFile(config['gene_annotation_gtf_file']['local_path']),
            ),
            kwargs={
                'overhang': config['star']['overhang'],
                'num_threads': config['star'].get('num_threads', 1),
            }
        )
    if 'tophat' in config:
        workflow.subworkflow(
            name='build_tophat_index',
            func=biowrappers.components.rna.tophat.workflow.create_tophat_transcriptome_index_workflow,
            args=(
                mgd.InputFile(config['ref_genome_fasta_file']['local_path']),
                mgd.InputFile(config['gene_annotation_gtf_file']['local_path']),
                mgd.OutputFile(config['tophat']['ref_genome_index']),
                mgd.OutputFile(config['tophat']['transcriptome_index']),
            ),
        )
    return workflow
