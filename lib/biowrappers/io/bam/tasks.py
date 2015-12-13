'''
Created on Nov 23, 2015

@author: Andrew Roth
'''
import pypeliner
import shutil

def convert_sam_to_bam(in_file, out_file):
    pass

def index_bam(bam_file, index_file):
    pypeliner.commandline.execute('samtools', 'index', bam_file)
    
    shutil.move(bam_file + '.bai', index_file)

def shuffle_bam(in_file, out_file):
    out_prefix = out_file.replace('.tmp', '').replace('.bam', '')
    
    pypeliner.commandline.execute('samtools', 'bamshuf', '-l', 9, in_file, out_prefix)
    
    shutil.move(out_prefix + '.bam', out_file)
    
def sort_bam(in_file, out_file, num_threads=1):
    out_prefix = out_file.replace('.tmp', '').replace('.bam', '')
    
    pypeliner.commandline.execute(
        'samtools', 
        'sort', 
        '-O', 'bam', 
        '-o', out_file,
        '-T', out_prefix,
        '-l', 9, 
        '-@', num_threads, 
        in_file)
