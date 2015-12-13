import pypeliner

import biowrappers.variant_calling.tri_nucleotide_context as tc

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
    
    pyp = pypeliner.app.Pypeline([tc.tasks], config)
    
    scheduler = pyp.sch
    
    tc.vcf_tric_nucleotide_annotation_pipeline(
        scheduler,
        args.ref_genome_fasta_file,
        args.vcf_file,
        args.out_file,
        split_size=args.split_size
    )
    
    pyp.run()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--ref_genome_fasta_file', required=True)
    
    parser.add_argument('--vcf_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--log_dir', default='./')
    
    parser.add_argument('--split_size', default=int(1e6), type=int)
        
    args = parser.parse_args()
    
    main(args)
    