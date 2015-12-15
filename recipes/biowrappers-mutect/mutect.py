#!/usr/bin/env python
import os
import subprocess
import sys

def main(args, extra_args):
    install_dir = os.path.dirname(os.path.realpath(__file__))

    java_dir = os.path.join(install_dir, 'java')
    
    if 'LD_LIBRARY_PATH' in os.environ:
        old_ld_path = os.environ['LD_LIBRARY_PATH']
    
    else:
        old_ld_path = ''
    
    # Make the sure the java VM can find the required libs
    os.environ['LD_LIBRARY_PATH'] = ':'.join((
        os.path.join(java_dir, 'lib', 'amd64', 'jli'),
        os.path.join(java_dir, 'jre', 'lib', 'amd64'),
        old_ld_path
    ))
    
    java_exe = os.path.join(java_dir, 'bin', 'java')

    # Fix memory issues
    os.environ['MALLOC_ARENA_MAX'] = '4'    
    
    cmd = [
        java_exe,
        '-Xms{0}g'.format(args.memory),
        '-Xmx{0}g'.format(args.memory),
        '-XX:PermSize=64M',
        '-XX:MaxPermSize=64M',
        '-XX:+UseConcMarkSweepGC',
        '-jar',
        os.path.join(install_dir, 'mutect.jar'),
        '--analysis_type', 'MuTect',
        '--input_file:normal', args.normal_bam_file,
        '--input_file:tumor', args.tumour_bam_file,
        '--reference_sequence', args.ref_genome_fasta_file,
        '--cosmic', args.cosmic_vcf_file,
        '--vcf', args.out_file,
        '--normal_sample_name', 'NORMAL',
        '--tumor_sample_name', 'TUMOUR',
        '-nt', 1,
        '-nct', 1
    ]
    
    if args.dbsnp_vcf_file is not None:
        cmd.extend(['--dbsnp', args.dbsnp_vcf_file])
        
    if args.region is not None:
        cmd.extend(['--intervals', args.region])
    
    cmd.extend(extra_args)
    
    run_cmd(cmd)

def run_cmd(cmd):
    cmd = [str(x) for x in cmd]

    process = subprocess.Popen(cmd)
    
    process.communicate()
    
    sys.exit(process.returncode)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='MuTect somatic variant caller')
#     
    parser.add_argument('--normal_bam_file', required=True)
     
    parser.add_argument('--tumour_bam_file', required=True)
     
    parser.add_argument('--ref_genome_fasta_file', required=True)
     
    parser.add_argument('--cosmic_vcf_file', required=True)
     
    parser.add_argument('--out_file', required=True)
     
    parser.add_argument('--dbsnp_vcf_file', default=None)
     
    parser.add_argument('--region', default=None)

    parser.add_argument('--memory', type=int, default=2)

    args, unknown = parser.parse_known_args()
    
    main(args, unknown)
