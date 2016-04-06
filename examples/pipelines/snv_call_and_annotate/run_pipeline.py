import os
import pypeliner
import yaml

from biowrappers.components.utils import make_directory
from biowrappers.components.variant_calling.utils import default_chromosomes
from biowrappers.pipelines.snv_call_and_annotate import call_and_annotate_pipeline 

import biowrappers.cli as cli

def main(args):
    make_directory(args.out_dir)
    
    with open(args.config_file) as fh:
        # Replace {ref_path_db} in config with desired path
        config_str = fh.read()
        
        config_str = config_str.format(ref_db_dir=args.ref_db_dir)
        
        # Load config
        config = yaml.load(config_str)
    
    if args.exome:
        config['strelka']['kwargs']['use_depth_thresholds'] = False
    
    tumour_bam_files = dict(zip(args.tumour_sample_ids, args.tumour_bam_files))
    
    raw_data_dir = os.path.join(args.out_dir, 'raw_data')
    
    results_dir = os.path.join(args.out_dir, 'results.h5')
    
    workflow = call_and_annotate_pipeline(
        config, 
        args.normal_bam_file, 
        tumour_bam_files, 
        raw_data_dir,
        results_dir, 
        chromosomes=args.chromosomes
    )
    
    pyp_config = cli.load_pypeliner_config(args)
            
    pyp = pypeliner.app.Pypeline([], pyp_config)
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--config_file', required=True)
    
    parser.add_argument('--ref_db_dir', required=True)
    
    parser.add_argument('--normal_bam_file', required=True)
    
    parser.add_argument('--tumour_bam_files', nargs='+', required=True)
    
    parser.add_argument('--tumour_sample_ids', nargs='+', required=True)
  
    parser.add_argument('--out_dir', required=True)
    
    parser.add_argument('--chromosomes', nargs='+', default=default_chromosomes)
    
    parser.add_argument('--exome', default=False, action='store_true')
    
    cli.add_pypeliner_args(parser)
    
    args = parser.parse_args()
    
    main(args)
