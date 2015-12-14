'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
import pypeliner

def run_classify(
    normal_bam_file,
    tumour_bam_file,
    ref_genome_fasta_file,
    region,
    out_file,
    min_normal_depth=1,
    min_tumour_depth=1,
    min_somatic_probability=0):
    
    cmd = [
        'bw-museq',
        'classify',
        '--normal_bam_file', normal_bam_file,
        '--tumour_bam_file', tumour_bam_file,
        '--ref_genome_fasta_file', ref_genome_fasta_file,
        '--out_file', out_file,
        '--min_normal_depth', min_normal_depth,
        '--min_tumour_depth', min_tumour_depth,
        '--min_somatic_probability', min_somatic_probability,
        '--region', region
    ]
    
    pypeliner.commandline.execute(*cmd)

def write_vcf(in_file, out_file, indel_threshold=0.05):
    
    cmd = [
        'bw-museq',
        'write_vcf',
        '--in_file', in_file,
        '--out_file', out_file,
        '--indel_threshold', indel_threshold
    ]
    
    pypeliner.commandline.execute(*cmd)
    