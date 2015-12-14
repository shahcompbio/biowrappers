import pypeliner

import biowrappers.cli as cli
import biowrappers.components.variant_calling.strelka as strelka

def main(args):
    config = cli.load_pypeliner_config(args)

    pyp = pypeliner.app.Pypeline([], config)

    indel_vcf_file = args.out_prefix + '.indel.vcf.gz'
    
    snv_vcf_file = args.out_prefix + '.snv.vcf.gz'
 
    workflow = strelka.strelka_pipeline(
        args.normal_bam_file, args.install_dir,
        args.tumour_bam_file, 
        args.ref_genome_fasta_file, 
        indel_vcf_file,
        snv_vcf_file,
        chromosomes=args.chromosomes,
        split_size=args.split_size,
        use_depth_thresholds=not args.no_depth_thresholds
    )
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    cli.add_normal_tumour_bam_variant_calling_args(parser)
    
    parser.add_argument('--out_prefix', required=True)
    
    cli.add_variant_calling_region_args(parser)

    parser.add_argument('--no_depth_thresholds', action='store_true', default=False)
    
    cli.add_pypeliner_args(parser)
        
    args = parser.parse_args()
    
    main(args)
    