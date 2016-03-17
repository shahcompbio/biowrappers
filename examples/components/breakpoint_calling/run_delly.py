import pypeliner

import biowrappers.cli as cli
import biowrappers.components.utils as utils
import biowrappers.components.breakpoint_calling.delly as delly

def main(args):
    utils.make_directory(args.raw_data_dir)

    config = cli.load_pypeliner_config(args)
    
    pyp = pypeliner.app.Pypeline(config=config)
    
    workflow = delly.delly_pipeline(
        args.normal_bam_file,
        cli.get_tumour_bam_file_dict(args),
        args.ref_genome_fasta_file,
        args.delly_excl_chrom,
        args.out_file,
        args.raw_data_dir,
    )
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()

    cli.add_normal_multiple_tumour_bam_args(parser)
    
    cli.add_ref_genome_arg(parser)

    parser.add_argument('--raw_data_dir', required=True)
        
    parser.add_argument('--out_file', required=True)
    
    cli.add_pypeliner_args(parser)
    
    parser.add_argument('--delly_excl_chrom', required=True)
        
    args = parser.parse_args()
    
    main(args)
