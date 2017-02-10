import os
import pypeliner
from pypeliner.workflow import Workflow

import biowrappers.components.io.hdf5.tasks as hdf5_tasks
import biowrappers.components.utils as utils
import biowrappers.components.io.download

import tasks


default_chromosomes = [str(a) for a in xrange(1, 23)] + ['X']


def create_ascat_workflow(
    seqdata_files,
    config,
    out_file,
    raw_data_dir,
    somatic_breakpoint_file=None,
    normal_id=None,
    **kwargs
):
    if normal_id is None:
        raise ValueError('ASCAT requires normal sample')

    normal_seqdata_file = seqdata_files[normal_id]
    tumour_seqdata_files = seqdata_files.copy()
    del tumour_seqdata_files[normal_id]

    results_files = os.path.join(raw_data_dir, 'results', 'sample_{sample_id}.h5')
    utils.make_parent_directory(results_files)

    workflow = Workflow()

    workflow.setobj(
        obj=pypeliner.managed.OutputChunks('sample_id'),
        value=tumour_seqdata_files.keys(),
    )

    workflow.transform(
        name='prepare_normal_data',
        ctx={'mem': 16, 'num_retry' : 3, 'mem_retry_increment' : 4},
        func=tasks.prepare_normal_data,
        args=(
            pypeliner.managed.InputFile(normal_seqdata_file),
            pypeliner.managed.TempOutputFile('Germline_LogR.txt'),
            pypeliner.managed.TempOutputFile('Germline_BAF.txt'),
            config,
        ),
    )

    workflow.transform(
        name='prepare_tumour_data',
        axes=('sample_id',),
        ctx={'mem': 20},
        func=tasks.prepare_tumour_data,
        args=(
            pypeliner.managed.InputFile('tumour_seqdata', 'sample_id', fnames=tumour_seqdata_files),
            pypeliner.managed.TempOutputFile('Germline_LogR.txt', 'sample_id'),
            pypeliner.managed.TempOutputFile('Germline_BAF.txt', 'sample_id'),
            config,
        ),
    )

    return workflow


