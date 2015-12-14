'''
Created on Nov 1, 2015

@author: Andrew Roth
'''

import pypeliner
import vcf

def run_vardict(
    normal_bam_file,
    tumour_bam_file,
    ref_genome_fasta_file,
    out_file,
    memory=8,
    min_allele_frequency=0.01,
    region=None,
    tumour_sample_name=None):
    
    cmd = [
        'bw-vardict',
        'run_vardict',
        '--normal_bam_file', normal_bam_file,
        '--tumour_bam_file', tumour_bam_file,
        '--ref_genome_fasta_file', ref_genome_fasta_file,
        '--out_file', out_file,
        '--memory', int(memory),
        '--min_allele_frequency', min_allele_frequency
    ]

    if region is not None:
        cmd.extend(['--region', region])
    
    if tumour_sample_name is not None:
        cmd.extend(['--tumour_sample_name', tumour_sample_name])
    
    pypeliner.commandline.execute(*cmd)
    
def run_vardict_test_somatic(in_file, out_file):

    cmd = [
        'bw-vardict',
        'test_somatic',
        '--in_file', in_file,
        '--out_file', out_file,
    ]
    
    pypeliner.commandline.execute(*cmd)
    
def run_vardict_var_to_vcf(
    in_file,
    out_file,
    min_allele_frequency=0.01):
    
    cmd = [
        'bw-vardict',
        'convert_to_vcf',
        '--in_file', in_file,
        '--out_file', out_file,
        '--min_allele_frequency', min_allele_frequency
    ]
    
    pypeliner.commandline.execute(*cmd)
    
def filter_vcf(in_file, out_file, variant_type):
    
    reader = vcf.Reader(filename=in_file)
    
    with open(out_file, 'w') as out_fh:
        writer = vcf.Writer(out_fh, reader)
        
        for record in reader:
            if record.INFO['STATUS'] != 'StrongSomatic':
                continue
            
            if len(record.FILTER) > 0:
                continue
            
            if (variant_type == 'indel') and (record.is_indel):
                writer.write_record(record)
            
            elif (variant_type == 'snv') and (record.is_snp):
                writer.write_record(record)                
