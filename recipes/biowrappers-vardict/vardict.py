#!/usr/bin/env python

import argparse
import os
import subprocess
import sys

install_dir = os.path.dirname(os.path.realpath(__file__))

bin_dir = os.path.join(install_dir, 'bin')

lib_dir = os.path.join(install_dir, 'lib')

scripts_dir = os.path.join(install_dir, 'scripts')

def cli():

    parser = argparse.ArgumentParser(description='VarDict variant caller')
    
    subparsers = parser.add_subparsers()
    
    run_vardict_parser = subparsers.add_parser('run_vardict')
    
    run_vardict_parser.add_argument('--normal_bam_file', required=True)
    
    run_vardict_parser.add_argument('--tumour_bam_file', required=True)
    
    run_vardict_parser.add_argument('--ref_genome_fasta_file', required=True)
    
    run_vardict_parser.add_argument('--out_file', required=True)
    
    run_vardict_parser.add_argument('--min_allele_frequency', default=0.01)
    
    run_vardict_parser.add_argument('--region', default=None)
    
    run_vardict_parser.add_argument('--tumour_sample_name', default=None)
    
    run_vardict_parser.add_argument('--memory', type=int, default=8)
    
    run_vardict_parser.set_defaults(func=run_vardict)
    
#---------------------------------------------------------------------------------------------------------------------- 
    
    test_somatic_parser = subparsers.add_parser('test_somatic')
    
    test_somatic_parser.add_argument('--in_file', required=True)
    
    test_somatic_parser.add_argument('--out_file', required=True)
    
    test_somatic_parser.set_defaults(func=test_somatic)
    
#---------------------------------------------------------------------------------------------------------------------- 
    
    convert_to_vcf_parser = subparsers.add_parser('convert_to_vcf')
    
    convert_to_vcf_parser.add_argument('--in_file', required=True)
    
    convert_to_vcf_parser.add_argument('--out_file', required=True)
    
    convert_to_vcf_parser.add_argument('--min_allele_frequency', default=0.01)
    
    convert_to_vcf_parser.set_defaults(func=convert_to_vcf)

#---------------------------------------------------------------------------------------------------------------------- 
    
    args = parser.parse_args()
    
    args.func(args)
    
def run_vardict(args):
    exe = os.path.join(bin_dir, 'VarDict')
    
    libs = ':'.join((
        os.path.join(lib_dir, 'VarDict-1.4.3.jar'),
        os.path.join(lib_dir, 'commons-cli-1.2.jar'),
        os.path.join(lib_dir, 'jregex-1.2_01.jar'),
        os.path.join(lib_dir, 'htsjdk-1.140.jar')
    ))
    
    os.environ['MALLOC_ARENA_MAX'] = '4'
    
    cmd = [
        'java',
        '-Xms{0}g'.format(args.memory),
        '-Xmx{0}g'.format(args.memory),
        '-XX:CompressedClassSpaceSize=64m',
        '-XX:MaxMetaspaceSize=128m',
        '-XX:+UseConcMarkSweepGC',
        '-classpath', libs,
        'com.astrazeneca.vardict.Main',
        '-G', args.ref_genome_fasta_file,
        '-f', args.min_allele_frequency,
        '-b', '"{0}|{1}"'.format(args.tumour_bam_file, args.normal_bam_file)
    ]
    
    if args.region is not None:
        cmd.extend(['-R', args.region])
    
    if args.tumour_sample_name is not None:
        cmd.extend(['-N', args.tumour_sample_name])
    
    run_cmd(cmd, args.out_file)

def test_somatic(args):
    exe = os.path.join(scripts_dir, 'testsomatic.R')
    
    cmd = [
        exe,
    ]
    
    run_cmd(cmd, args.out_file, in_file=args.in_file)

def convert_to_vcf(args):
    exe = os.path.join(scripts_dir, 'var2vcf_somatic.pl')
    
    cmd = [
        exe,
        '-f', args.min_allele_frequency,
    ]
    
    run_cmd(cmd, args.out_file, in_file=args.in_file)

def run_cmd(cmd, out_file, in_file=None):
    cmd = [str(x) for x in cmd]
    
    with open(out_file, 'w') as out_fh:
        if in_file is not None:
            with open(in_file, 'r') as in_fh:
                process = subprocess.Popen(cmd, stdin=in_fh, stdout=out_fh, shell=False)
           
            
        else:           
            process = subprocess.Popen(cmd, stdout=out_fh, shell=False)
        
        process.communicate()
  
    sys.exit(process.returncode)
  
if __name__ == '__main__':
    cli()
    
