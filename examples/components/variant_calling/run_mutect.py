import pypeliner

import biowrappers.cli as cli
import biowrappers.components.variant_calling.mutect as mutect


def main(args):
    config = cli.load_pypeliner_config(args)

    pyp = pypeliner.app.Pypeline([], config)

    workflow = mutect.create_mutect_workflow(
        args.normal_bam_file,
        args.tumour_bam_file,
        args.ref_genome_fasta_file,
        args.cosmic_file,
        args.dbsnp_file,
        args.out_file,
        chromosomes=args.chromosomes,
        split_size=args.split_size
    )

    pyp.run(workflow)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    cli.add_normal_tumour_bam_args(parser)

    cli.add_ref_genome_arg(parser)

    parser.add_argument('--cosmic_file', required=True)

    parser.add_argument('--dbsnp_file', required=True)

    parser.add_argument('--out_file', required=True)

    cli.add_variant_calling_region_args(parser)

    cli.add_pypeliner_args(parser)

    args = parser.parse_args()

    main(args)
