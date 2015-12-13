import os
import pypeliner

import biowrappers.variant_calling.snpeff as snpeff

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
    
    pyp = pypeliner.app.Pypeline([snpeff.tasks], config)
    
    scheduler = pyp.sch
    
    snpeff.snpeff_pipeline(
        scheduler,
        args.install_dir,
        args.target_vcf_file,   
        args.out_file,
        data_base=args.data_base,
        split_size=args.split_size
    )
    
    pyp.run()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--install_dir', required=True)
    
    parser.add_argument('--target_vcf_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--data_base', default='GRCh37.75')
    
    parser.add_argument('--log_dir', default='./')
    
    parser.add_argument('--split_size', default=int(1e3), type=int)
    
    args = parser.parse_args()
    
    main(args)
    