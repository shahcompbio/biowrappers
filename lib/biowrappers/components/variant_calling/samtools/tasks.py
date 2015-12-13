'''
Created on Nov 2, 2015

@author: Andrew Roth
'''
import pypeliner

def run_samtools_variant_calling(
        bam_file,
        ref_genome_fasta_file,
        out_file,
        max_depth=int(1e7),
        min_bqual=0,
        min_depth=0,
        min_mqual=0,
        region=None):
    
    mpileup_cmd = [
        'samtools',
        'mpileup',
        '-u',
        '-f', ref_genome_fasta_file,
        '-Q', min_bqual,
        '-q', min_mqual,
        bam_file
    ] 
    
    if region is not None:
        mpileup_cmd.extend(['-r', region])
    
    bcf_cmd_1 = [
        'bcftools',
        'view',
        '-bvcg',
        '-'  
    ]
    
    bcf_cmd_2 = [
        'bcftools',
        'view',
        '-'
    ]
    
    vcf_cmd = [
        'vcfutils.pl',
        'varFilter',
        '-d', min_depth,
        '-D', max_depth
    ]
    
    cmd = []
    
    cmd.extend(mpileup_cmd)
    cmd.append('|')
    cmd.extend(bcf_cmd_1)
    cmd.append('|')
    cmd.extend(bcf_cmd_2)
    cmd.append('|')
    cmd.extend(vcf_cmd)
    cmd.extend(['>', out_file])
    
    pypeliner.commandline.execute(*cmd)
