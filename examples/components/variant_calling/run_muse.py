import pypeliner

import biowrappers.variant_calling.muse as muse

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
    
    pyp = pypeliner.app.Pypeline([muse.tasks], config)
    
    scheduler = pyp.sch
    
    muse.muse_pipeline(
        scheduler, 
        args.normal_bam_file, 
        args.tumour_bam_file, 
        args.ref_genome_fasta_file, 
        args.out_file,
        data_type=args.data_type,
        chromosomes=['21',],
        split_size=args.split_size
    )
    
    pyp.run()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--normal_bam_file', required=True)
    
    parser.add_argument('--tumour_bam_file', required=True)
    
    parser.add_argument('--ref_genome_fasta_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--data_type', choices=['exome', 'wgsss'], default='wgss')
    
    parser.add_argument('--log_dir', default='./')
    
    parser.add_argument('--split_size', default=int(1e6), type=int)
        
    args = parser.parse_args()
    
    main(args)
    