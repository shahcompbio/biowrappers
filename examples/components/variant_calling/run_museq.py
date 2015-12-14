import pypeliner

import biowrappers.cli as cli
import biowrappers.components.variant_calling.museq as museq

def main(args):
    config = cli.load_pypeliner_config(args)
    
    pyp = pypeliner.app.Pypeline([], config)
    
    workflow = museq.museq_pipeline(
        args.model_file,
        args.normal_bam_file,
        args.tumour_bam_file,
        args.ref_genome_fasta_file,
        args.out_file,
        chromosomes=args.chromosomes,
        split_size=args.split_size
    )
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--model_file', required=True)
    
    cli.add_normal_tumour_bam_variant_calling_args(parser)
    
    parser.add_argument('--out_file', required=True)
    
    cli.add_variant_calling_region_args(parser)
    
    cli.add_pypeliner_args(parser)
    
    args = parser.parse_args()
    
    main(args)
