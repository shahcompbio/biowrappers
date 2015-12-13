import pypeliner

from biowrappers.variant_calling.utils import default_chromosomes

import biowrappers.variant_calling.snv_allele_counts as snv_allele_counts

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
    
    tumour_bam_files = dict(zip(args.tumour_samples, args.tumour_bam_files))
    
    pyp = pypeliner.app.Pypeline([snv_allele_counts.tasks], config)
    
    scheduler = pyp.sch
    
    snv_allele_counts.snv_variant_position_counts_pipeline(
        scheduler,
        args.normal_bam_file,
        tumour_bam_files,
        args.out_file,
        chromosomes=args.chromosomes,
        count_duplicates=args.count_duplicates,
        min_bqual=args.min_bqual,
        min_mqual=args.min_mqual,
        min_normal_depth=args.min_normal_depth,
        min_tumour_depth=args.min_tumour_depth,
        min_variant_depth=args.min_variant_depth,
        report_strand_counts=args.report_strand_counts,
        split_size=args.split_size)
    
    pyp.run()

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--normal_bam_file', required=True)
    
    parser.add_argument('--tumour_bam_files', nargs='+', required=True)
    
    parser.add_argument('--tumour_samples', nargs='+', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--count_duplicates', action='store_true', default=False)
    
    parser.add_argument('--chromosomes', default=default_chromosomes, nargs='+')
    
    parser.add_argument('--log_dir', default='./')
    
    parser.add_argument('--min_bqual', default=30, type=int)
    
    parser.add_argument('--min_mqual', default=30, type=int)
    
    parser.add_argument('--min_normal_depth', default=0, type=int)
    
    parser.add_argument('--min_tumour_depth', default=0, type=int)
    
    parser.add_argument('--min_variant_depth', default=0, type=int)
    
    parser.add_argument('--report_strand_counts', action='store_true', default=False)
    
    parser.add_argument('--split_size', default=int(1e6), type=int)
    
    args = parser.parse_args()
    
    main(args)
    
