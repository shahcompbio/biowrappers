import pypeliner

import biowrappers.cli as cli
import biowrappers.components.variant_calling.museq as museq

def main(args):
    config = cli.load_pypeliner_config(args)
    
    pyp = pypeliner.app.Pypeline([], config)
    
    workflow = museq.museq_pipeline(
        args.normal_bam_file,
        args.tumour_bam_files,
        args.ref_genome_fasta_file,
        args.out_file,
        chromosomes=args.chromosomes,
        indel_threshold=args.indel_threshold,
        chunk_size=args.chunk_size,
        min_normal_depth=args.min_normal_depth,
        min_tumour_depth=args.min_tumour_depth,
        min_somatic_probability=args.min_somatic_probability,
        split_size=args.split_size
    )
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()

    cli.add_normal_multiple_tumour_bam_variant_calling_args(parser)
    
    parser.add_argument('--out_file', required=True)
    
    cli.add_variant_calling_region_args(parser)
    
    cli.add_pypeliner_args(parser)
    
    parser.add_argument('--chunk_size', default=int(1e5), type=int)
    
    parser.add_argument('--indel_threshold', default=0.05, type=float)
    
    parser.add_argument('--min_normal_depth', default=1, type=int)
    
    parser.add_argument('--min_tumour_depth', default=1, type=int)
    
    parser.add_argument('--min_somatic_probability', default=0.5, type=float)
    
    args = parser.parse_args()
    
    main(args)
