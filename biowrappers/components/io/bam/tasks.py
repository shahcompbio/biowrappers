'''
Created on Feb 3, 2016

@author: Andrew Roth
'''
import glob
import gzip
import pypeliner
import os
import shutil
import tempfile
import time

from biowrappers.components.utils import flatten_input

def collate(in_file, out_file):
    
    out_prefix = out_file + '.shuffle'
    
    pypeliner.commandline.execute('samtools', 'collate', in_file, out_prefix)
    
    time.sleep(1)
    
    shutil.move(out_prefix + '.bam', out_file)

def convert_to_fastqs(in_file, read_files, tmp_dir, split_size=int(1e7)):
    
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    
    # Clean tmp area
    for file_name in glob.glob(os.path.join(tmp_dir, '*')):
        os.unlink(file_name)
    
    prefix = {}
    
    read_ids = [1, 2]
    
    for r_id in read_ids:
        prefix[r_id] = os.path.join(tmp_dir, 'read_{0}'.format(r_id))
    
    pypeliner.commandline.execute(
        'bamtofastq',
        'F={0}'.format(prefix[1]),
        'F2={0}'.format(prefix[2]),
        'gz=1',
        'filename={0}'.format(in_file),
        'split={0}'.format(split_size),
        '>', '/dev/null'
    )
    
    num_files = len(glob.glob(prefix[1] + '*.gz'))
    
    last_file_idx = num_files - 1
    
    num_fn_digits = 6
    
    len_last_file_idx = len(str(last_file_idx))
    
    first_file_idx = '0' * num_fn_digits
    
    last_file_idx = '0' * (num_fn_digits - len_last_file_idx) + str(last_file_idx)

    # Merge last file into first so all files have at least minimum number of reads
    if num_files > 1:
        for r_id in read_ids:
            first_file = os.path.join(tmp_dir, 'read_{0}_{1}.gz'.format(r_id, first_file_idx))
             
            last_file = os.path.join(tmp_dir, 'read_{0}_{1}.gz'.format(r_id, last_file_idx))
             
            with gzip.GzipFile(last_file, 'r') as in_fh, gzip.GzipFile(first_file, 'a') as out_fh:
                for line in in_fh:
                    out_fh.write(line)
             
            os.unlink(last_file)
             
    # Move files over
    for r_id in read_ids:
        for tmp_file in glob.glob(prefix[r_id] + '*.gz'):
            idx = int(os.path.splitext(os.path.basename(tmp_file))[0].split('_')[-1])
             
            out_file = read_files[r_id][idx]
             
            shutil.move(tmp_file, out_file)

def mark_duplicates(
    in_files, 
    out_file, 
    compression_level=9,
    hash_table_size=262144,    
    io_buffer_size=128,
    num_threads=1, 
    overflow_list_size=200000,
    tmp_dir=None):
    
    try:
        if tmp_dir is None:
            tmp_dir = tempfile.mkdtemp()
            
            clean_up = True
        
        else:
            clean_up = False
        
        cmd = [
            'sambamba',
            'markdup',
            '-p',
            '-l', compression_level,
            '-t', num_threads,
            '--io-buffer-size={0}'.format(io_buffer_size),
            '--hash-table-size={0}'.format(hash_table_size),
            '--overflow-list-size={0}'.format(overflow_list_size),
            '--tmpdir', tmp_dir,
        ]
        
        cmd.extend(flatten_input(in_files))
        
        cmd.append(out_file)
        
        pypeliner.commandline.execute(*cmd)
    
    finally:
        if clean_up:
            shutil.rmtree(tmp_dir)

def merge(
    in_files,
    out_file,
    attach_read_group_from_file_name=False,
    header_file=None,
    compression_level=9,
    num_compression_threads=0):
    
    cmd = [
        'samtools',
        'merge',
        '-c',
        '-p',
        '-f',
        '-l', compression_level,
        '-@', num_compression_threads,
    ]
    
    if attach_read_group_from_file_name:
        cmd.append('-r')
        
    if header_file is not None:
        cmd.extend(['-h', header_file])
    
    cmd.append(out_file)
    
    for file_name in flatten_input(in_files):
        cmd.append(file_name)
        
    pypeliner.commandline.execute(*cmd)
    
def sort(in_file, out_file, max_mem='2G', name_sort=False, compression_level=9, num_threads=1):

    cmd = [
        'samtools',
        'sort',
        '-l', compression_level,
        '-m', max_mem,
        '-o', out_file,
        '-@', num_threads,
    ]
     
    if name_sort:
        cmd.append('-n')
     
    cmd.append(in_file)
     
    pypeliner.commandline.execute(*cmd)

