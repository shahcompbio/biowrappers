'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
import pypeliner
import shutil

def run_muse_call(
    normal_bam_file,
    tumour_bam_file,
    ref_genome_fasta_file,
    out_file,
    region=None):
    
    tmp_file = '.'.join(out_file.split('.')[:-3])

    cmd = [
        'muse',
        'call',
        '-f', ref_genome_fasta_file,
        '-O', tmp_file,
        tumour_bam_file,
        normal_bam_file
    ]
    
    if region is not None:
        cmd.extend(['-r', region, ])
    
    pypeliner.commandline.execute(*cmd)
    
    shutil.move(tmp_file + '.MuSE.txt', out_file)
    
def run_muse_sump(
    in_file,
    out_file,
    exome=False):

    cmd = [
        'muse',
        'sump',
        '-I', in_file,
        '-O', out_file
    ]
    
    if exome:
        cmd.append('-E')
    
    else:
        cmd.append('-G')
    
    pypeliner.commandline.execute(*cmd)           
