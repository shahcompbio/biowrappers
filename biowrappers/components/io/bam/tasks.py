'''
Created on Feb 3, 2016

@author: Andrew Roth
'''
import gzip
import glob
import pypeliner
import os
import shutil
import time

def collate_bam(in_file, out_file):
    
    out_prefix = out_file + '.shuffle'
    
    pypeliner.commandline.execute('samtools', 'collate', in_file, out_prefix)
    
    time.sleep(1)
    
    shutil.move(out_prefix + '.bam', out_file)

def split_bam_to_fastq(in_file, read_files, tmp_dir, split_size=int(1e7)):
    print in_file, tmp_dir
    
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
        '>/dev/null'
    )
    
    num_files = len(glob.glob(prefix[1] + '*.gz'))
        
    last_file_idx = num_files - 1
    
    num_fn_digits = 6
    
    len_last_file_idx = len(str(last_file_idx))
    
    first_file_idx = '0' * num_fn_digits
    
    last_file_idx = '0' * (num_fn_digits - len_last_file_idx) + str(last_file_idx)

    # Merge last file into first so all files have at least minimum number of reads
    if last_file_idx != 0:
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
