import os
import yaml
import pypeliner
from pypeliner.workflow import Workflow

import biowrappers.components.utils
import biowrappers.components.io.download
import biowrappers.pipelines.realignment
import biowrappers.components.io.bam.tasks
import biowrappers.pipelines.breakpoint_call_and_annotate
import biowrappers.pipelines.copy_number

import tasks


def main(args):
    biowrappers.components.utils.make_directory(args.out_dir)
    
    with open(args.config_file) as config_file:
        config_text = config_file.read()
    config_text = config_text.format(out_dir=args.out_dir, ref_db_dir=args.ref_db_dir)
    config = yaml.load(config_text)

    pypeliner_args = vars(args)
    pypeliner_args['tmpdir'] = os.path.join(args.out_dir, 'pipeline')

    pyp = pypeliner.app.Pypeline(modules=[tasks], config=pypeliner_args)

    download_urls = {}

    for sample in ('tumour', 'normal'):
        lanes = config['lanes'][sample]

        for lane in lanes:
            download_urls[(sample, lane)] = config['lanes'][sample][lane]['url']

    raw_lane_template = os.path.join(args.out_dir, 'lanes', 'raw', '{lane}.bam')

    realigned_lane_template = os.path.join(args.out_dir, 'lanes', 'realigned', '{lane}.bam')
    sample_bam_template = os.path.join(args.out_dir, 'lanes', '{sample}.bam')

    workflow = Workflow(default_ctx={'mem': 8})

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('url', 'sample', 'lane'),
        value=download_urls,
    )

    workflow.subworkflow(
        name='download_lanes',
        axes=('sample', 'lane'),
        func=biowrappers.components.io.download.create_download_workflow,
        args=(
            pypeliner.managed.TempInputObj('url', 'sample', 'lane'),
            pypeliner.managed.OutputFile('raw_lane', 'sample', 'lane', template=raw_lane_template),
        )
    )

    workflow.subworkflow(
        name='realign_lanes',
        axes=('sample', 'lane'),
        func=biowrappers.pipelines.realignment.realignment_pipeline,
        args=(
            config['realignment'],
            pypeliner.managed.InputFile('raw_lane', 'sample', 'lane', template=raw_lane_template),
            pypeliner.managed.OutputFile('realigned_lane', 'sample', 'lane', template=realigned_lane_template),
        )
    )

    workflow.transform(
        name='merge_and_markdups',
        axes=('sample',),
        func=biowrappers.components.io.bam.tasks.mark_duplicates,
        args=(
            pypeliner.managed.InputFile('realigned_lane', 'sample', 'lane', template=realigned_lane_template),
            pypeliner.managed.OutputFile('bam', 'sample', template=sample_bam_template),
        ),
        kwargs={
            'tmp_dir': pypeliner.managed.TempSpace('markdup_temp', 'sample')
        }
    )

    normal_bam_file = sample_bam_template.format(sample='normal')
    tumour_bam_file = sample_bam_template.format(sample='tumour')

    breakpoint_raw_data_dir = os.path.join(args.out_dir, 'breakpoints', 'raw')
    breakpoint_results_file = os.path.join(args.out_dir, 'breakpoints', 'results.h5')

    workflow.subworkflow(
        name='breakpoint_call_and_annotate',
        func=biowrappers.pipelines.breakpoint_call_and_annotate.call_and_annotate_pipeline,
        args=(
            config,
            pypeliner.managed.InputFile(normal_bam_file),
            {'tumour': pypeliner.managed.InputFile(tumour_bam_file)},
            pypeliner.managed.Template(os.path.join(breakpoint_raw_data_dir)),
            pypeliner.managed.OutputFile(breakpoint_results_file),
        ),
    )

    somatic_breakpoints_file = os.path.join(args.out_dir, 'somatic_breakpoints.tsv')

    workflow.transform(
        name='extract_somatic_breakpoint',
        ctx={'mem':4},
        func=tasks.extract_somatic_breakpoint,
        args=(
            pypeliner.managed.InputFile(breakpoint_results_file),
            pypeliner.managed.OutputFile(somatic_breakpoints_file),
            config,
        )
    )

    copy_number_raw_data_dir = os.path.join(args.out_dir, 'copy_number', 'raw')
    breakpoint_results_file = os.path.join(args.out_dir, 'copy_number', 'results.h5')

    workflow.subworkflow(
        name='copy_number_call_and_annotate',
        func=biowrappers.pipelines.copy_number.call_and_annotate_pipeline,
        args=(
            config,
            pypeliner.managed.InputFile(normal_bam_file),
            {'tumour': pypeliner.managed.InputFile(tumour_bam_file)},
            pypeliner.managed.InputFile(somatic_breakpoints_file),
            copy_number_raw_data_dir,
            pypeliner.managed.OutputFile(breakpoint_results_file),
        ),
    )

    pyp.run(workflow)



if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()

    pypeliner.app.add_arguments(parser)

    parser.add_argument('--config_file', required=True)
        
    parser.add_argument('--ref_db_dir', required=True)

    parser.add_argument('--out_dir', required=True)

    args = parser.parse_args()
    
    main(args)
