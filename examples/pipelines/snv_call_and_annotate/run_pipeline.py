import pypeliner
import yaml

from biowrappers.components.utils import make_directory
from biowrappers.components.variant_calling.utils import default_chromosomes
from biowrappers.pipelines.snv_call_and_annotate import call_and_annotate_pipeline 

def main(args):
    make_directory(args.out_dir)
    
    with open(args.config_file) as fh:
        config = yaml.load(fh)
    
    if args.exome:
        config['strelka']['kwargs']['use_depth_thresholds'] = False
    
    tumour_bam_files = dict(zip(args.tumour_sample_ids, args.tumour_bam_files))
    
    workflow = call_and_annotate_pipeline(
        config, 
        args.normal_bam_file, 
        tumour_bam_files, 
        args.ref_genome_fasta_file, 
        args.out_dir, 
        chromosomes=args.chromosomes
    )
        
    if args.log_dir is not None:
        config['pypeliner']['tmpdir'] = args.log_dir
        
    pyp = pypeliner.app.Pypeline([], config['pypeliner'])
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--config_file', required=True)
    
    parser.add_argument('--normal_bam_file', required=True)
    
    parser.add_argument('--tumour_bam_files', nargs='+', required=True)
    
    parser.add_argument('--tumour_sample_ids', nargs='+', required=True)
    
    parser.add_argument('--ref_genome_fasta_file', required=True)
    
    parser.add_argument('--out_dir', required=True)
    
    parser.add_argument('--chromosomes', nargs='+', default=default_chromosomes)
    
    parser.add_argument('--exome', default=False, action='store_true')
    
    parser.add_argument('--log_dir', default='./')
    
    args = parser.parse_args()
    
    main(args)
