import pypeliner

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.variant_calling.mutect as mutect

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
    
    pyp = pypeliner.app.Pypeline([], config)
    
    workflow = mutect.mutect_pipeline(
        args.normal_bam_file, 
        args.tumour_bam_file, 
        args.ref_genome_fasta_file,
        args.cosmic_file,
        args.dbsnp_file,
        args.out_file, 
        chromosomes=args.chromosomes,
        memory=args.memory,
        split_size=args.split_size
    )
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--normal_bam_file', required=True)
    
    parser.add_argument('--tumour_bam_file', required=True)
    
    parser.add_argument('--ref_genome_fasta_file', required=True)
    
    parser.add_argument('--cosmic_file', required=True)
    
    parser.add_argument('--dbsnp_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--chromosomes', nargs='+', default=default_chromosomes)
    
    parser.add_argument('--log_dir', default='./')
    
    parser.add_argument('--memory', default=4, type=int)
    
    parser.add_argument('--split_size', default=int(1e6), type=int)
        
    args = parser.parse_args()
    
    main(args)
