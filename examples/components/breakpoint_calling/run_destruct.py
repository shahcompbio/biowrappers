import yaml
import pypeliner

import biowrappers.cli as cli
import biowrappers.components.utils as utils
import biowrappers.components.breakpoint_calling.destruct as destruct


def main(args):
    utils.make_directory(args.raw_data_dir)

    destruct_config = {}
    if args.destruct_config is not None:
        with open(args.destruct_config) as fh:
            destruct_config = yaml.load(fh)

    pypeliner_config = cli.load_pypeliner_config(args)

    pyp = pypeliner.app.Pypeline(config=pypeliner_config)

    workflow = destruct.destruct_pipeline(
        args.normal_bam_file,
        cli.get_tumour_bam_file_dict(args),
        destruct_config,
        args.ref_data_dir,
        args.out_file,
        args.raw_data_dir,
    )

    pyp.run(workflow)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    cli.add_normal_multiple_tumour_bam_args(parser)
    cli.add_pypeliner_args(parser)

    parser.add_argument('--out_file', required=True)
    parser.add_argument('--raw_data_dir', required=True)
    parser.add_argument('--ref_data_dir', required=True)
    parser.add_argument('--destruct_config', default=None)

    args = parser.parse_args()

    main(args)
