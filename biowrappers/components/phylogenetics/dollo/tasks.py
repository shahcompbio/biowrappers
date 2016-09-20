'''
Created on Feb 29, 2016

@author: Andrew Roth
'''
import glob
import gzip
import pandas as pd
import pypeliner
import os
import shutil
import pickle

def create_trees(in_file, out_files, tmp_dir):
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        
    cmd = [
        'PyDollo',
        'create_trees',
        '--log_likelihoods_file', in_file,
        '--out_dir', tmp_dir
    ]
    
    pypeliner.commandline.execute(*cmd)
    
    # Move files from tmp space
    for tmp_file in glob.glob(os.path.join(tmp_dir, '*.pickle.gz')):
        idx = int(os.path.basename(tmp_file).split('.')[0])
         
        out_file = out_files[idx]
         
        shutil.move(tmp_file, out_file)
    
    os.rmdir(tmp_dir)

def compute_tree_log_likelihoods(
        log_likelihoods_file,
        tree_file, 
        out_file,
        max_probability_of_loss=0.5,
        min_probability_of_loss=0.0, 
        probability_of_loss=None):
    
    cmd = [
        'PyDollo',
        'compute_tree_log_likelihood',
        '--log_likelihoods_file', log_likelihoods_file,
        '--tree_file', tree_file,
        '--out_file', out_file,
    ]
     
    
    if probability_of_loss is None:
        cmd.extend([
            '--max_probability_of_loss', max_probability_of_loss,
            '--min_probability_of_loss', min_probability_of_loss,
        ])
    
    else:
        cmd.extend(['--probability_of_loss', probability_of_loss])
    
    pypeliner.commandline.execute(*cmd)

def select_ml_tree(optimal_files, out_file):
    tree_log_likelihoods = []
    
    for tree_id in optimal_files:
        results = pickle.load(gzip.open(optimal_files[tree_id], 'r'))
        tree_log_likelihoods.append((results['log_likelihood'], tree_id))
        
    ml_tree_id = sorted(tree_log_likelihoods)[-1][1]
    
    shutil.copyfile(optimal_files[ml_tree_id], out_file)
        
