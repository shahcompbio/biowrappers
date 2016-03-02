'''
Created on Dec 13, 2015

@author: Andrew Roth
'''
def add_pypeliner_args(parser):
    parser.add_argument('-c', '--cleanup_tmp_files', default=True, action='store_false')
    
    parser.add_argument('-l', '--log_dir', default='./')
    
    parser.add_argument('-m', '--max_jobs', default=1)
    
    parser.add_argument('-n', '--native_spec', default='-V -q all.q -l mem_token={mem}G,mem_free={mem}G,h_vmem={mem}G')
    
    parser.add_argument('-s', '--submit_method', choices=['asyncqsub', 'drmaa', 'local'], default='local')

    parser.add_argument('-rp', '--repopulate', action = "store_true")

def add_variant_calling_region_args(parser):
    from biowrappers.components.variant_calling.utils import default_chromosomes

    parser.add_argument('--chromosomes', nargs='+', default=default_chromosomes)
    
    parser.add_argument('--split_size', default=int(1e7), type=int)
    
def add_normal_tumour_bam_args(parser):
    parser.add_argument('-nb', '--normal_bam_file', required=True)
    
    parser.add_argument('-tb', '--tumour_bam_file', required=True)

def add_normal_multiple_tumour_bam_args(parser):
    parser.add_argument('-nb', '--normal_bam_file', required=True)
    
    parser.add_argument('-tb', '--tumour_bam_files', nargs='+', required=True)

def add_ref_genome_arg(parser):
    parser.add_argument('-rg', '--ref_genome_fasta_file', required=True)    

def load_pypeliner_config(args):
    config = {
        'tmpdir' : args.log_dir,
        'pretend' : False,
        'submit' : args.submit_method,
        'nativespec' : args.native_spec,
        'maxjobs' : args.max_jobs,
        'nocleanup' : not args.cleanup_tmp_files
    }
    
    return config
    
