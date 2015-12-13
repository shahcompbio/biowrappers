import pypeliner

import biowrappers.variant_calling.snv_allele_counts as snv_allele_counts

def main(args):
    native_spec = '-V -q all.q -l mem_token={mem}G,mem_free={mem}G,h_vmem={mem}G'
    
    config = {
        'tmpdir' : args.log_dir,
        'pretend' : False,
        'submit' : 'asyncqsub',
        'nativespec' : native_spec,
        'maxjobs' : 100,
        'nocleanup' : True
    }
    
    pyp = pypeliner.app.Pypeline([snv_allele_counts.tasks, snv_allele_counts.vcf_tasks], config)
    
    scheduler = pyp.sch
    
    snv_allele_counts.snv_allele_counts_for_vcf_targets_pipeline(
        scheduler, 
        args.bam_file, 
        args.vcf_file,  
        args.out_file,
        count_duplicates=args.count_duplicates,
        min_bqual=args.min_bqual,
        min_mqual=args.min_mqual,
        split_size=args.split_size)
    
    pyp.run()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--bam_file', required=True)
    
    parser.add_argument('--vcf_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--count_duplicates', action='store_true', default=False)
    
    parser.add_argument('--min_bqual', default=30, type=int)
    
    parser.add_argument('--min_mqual', default=30, type=int)
    
    parser.add_argument('--log_dir', default='./')
    
    parser.add_argument('--split_size', default=int(1e4), type=int)
    
    args = parser.parse_args()
    
    main(args)
    