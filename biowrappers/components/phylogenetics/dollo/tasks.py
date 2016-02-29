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

def build_results_file(optimal_files, out_file, grid_search_files=None):
    data = []
    
    for tree_id in optimal_files:
        df = pd.read_csv(optimal_files[tree_id], compression='gzip', sep='\t')

        row = {
               'tree_id' : tree_id,
               'probability_of_loss' : df.at[0, 'probability_of_loss'],
               'log_likelihood' : df.at[0, 'log_likelihood']
               }
        
        data.append(row)
    
    if grid_search_files is not None:
        for key in grid_search_files:
            df = pd.read_csv(grid_search_files[key], compression='gzip', sep='\t')
    
            row = {
                   'tree_id' : key[0],
                   'probability_of_loss' : df.at[0, 'probability_of_loss'],
                   'log_likelihood' : df.at[0, 'log_likelihood']
                   }
            
            data.append(row)
    
    data = pd.DataFrame(data, columns=['tree_id', 'probability_of_loss', 'log_likelihood'])
    
    data = data.sort('log_likelihood', ascending=False)
    
    with gzip.GzipFile(out_file, 'w') as out_fh:
        data.to_csv(out_fh, index=False, sep='\t')
        
def find_ml_tree(search_file, tree_files, out_file):
    print search_file
    
    df = pd.read_csv(search_file, compression='gzip', sep='\t')
    
    row = df.iloc[df.log_likelihood.idxmax()]
    
    best_file = tree_files[int(row['tree_id'])]
    
    shutil.copyfile(best_file, out_file)