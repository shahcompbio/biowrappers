'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
import os
import pypeliner
import time

def run_mutect(
    normal_bam_file,
    tumour_bam_file,
    ref_genome_fasta_file,
    cosmic_vcf_file,
    dbsnp_vcf_file,
    region,
    out_file):
    
    cmd = [
        'bw-mutect',
        '--normal_bam_file', normal_bam_file,
        '--tumour_bam_file', tumour_bam_file,
        '--ref_genome_fasta_file', ref_genome_fasta_file,
        '--cosmic_vcf_file', cosmic_vcf_file,
        '--dbsnp_vcf_file', dbsnp_vcf_file,
        '--out_file', out_file,
        '--region', region,
        '>', '/dev/null' # Avoid dumping mutect tsv file to stdout
    ]

    pypeliner.commandline.execute(*cmd)

    idx_file_name = out_file + '.idx'
    
    # Guard against slow file system
    time.sleep(1)
    
    os.unlink(idx_file_name)
