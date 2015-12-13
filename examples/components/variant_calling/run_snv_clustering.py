import os
import pypeliner

from biowrappers.variant_calling.utils import default_chromosomes

import biowrappers.variant_calling.snv_clustering as snv_clustering

def main(args):
    native_spec = '-V -q all.q -l mem_token={mem}G,mem_free={mem}G,h_vmem={mem}G'
    
    config = {
        'tmpdir' : args.log_dir,
        'pretend' : False,
        'submit' : 'asyncqsub',
        'nativespec' : native_spec,
        'maxjobs' : 100,
        'nocleanup' : False
    }
    
    pyp = pypeliner.app.Pypeline([snv_clustering.tasks], config)
    
    scheduler = pyp.sch
    
    snv_clustering.snv_allele_count_clustering_pipeline(
        scheduler, 
        args.counts_file,
        args.out_file,
        num_restarts=args.num_restarts)
    
    pyp.run()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--counts_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--log_dir', default='./')
    
    parser.add_argument('--num_restarts', default=10, type=int)
    
    args = parser.parse_args()
    
    main(args)
    