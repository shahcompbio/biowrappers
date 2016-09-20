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
import collections

import dollo.trees
import dollo.run


def create_tree_groups(log_likelihoods_file, trees_files, max_tree_groups=1000):
    ll = pd.read_csv(log_likelihoods_file, compression='gzip', sep='\t')
    
    leaf_names = ll['site_id'].unique()
    
    tree_groups = collections.defaultdict(dict)
    for tree_id, tree in enumerate(dollo.trees.enumerate_labeled_trees(leaf_names)):
        tree_groups[tree_id % max_tree_groups][tree_id] = tree

    for idx, trees in tree_groups.iteritems():
        with gzip.open(trees_files[idx], 'w') as f:
            pickle.dump(trees, f)


def compute_tree_log_likelihoods(
        log_likelihoods_file,
        trees_file, 
        results_file,
        max_probability_of_loss=0.5,
        min_probability_of_loss=0.0, 
        probability_of_loss=None):
    
    ll = pd.read_csv(log_likelihoods_file, compression='gzip', sep='\t')
    trees = pickle.load(gzip.open(trees_file))
    
    args = {
        'min_probability_of_loss': min_probability_of_loss,
        'max_probability_of_loss': max_probability_of_loss,
        'probability_of_loss': probability_of_loss,
    }
    
    results = {}
    for tree_id, tree in trees.iteritems():
        result = dollo.run.compute_tree_log_likelihood(ll, tree, **args)
        results[tree_id] = result
        
    with gzip.open(results_file, 'w') as f:
        pickle.dump(results, f)


def select_ml_tree(results_files, ml_result_file, search_file):
    trees = {}
    tree_log_likelihoods = []
    
    for tree_group_id in results_files:
        results = pickle.load(gzip.open(results_files[tree_group_id], 'r'))
        
        for tree_id, result in results.iteritems():
            tree_log_likelihoods.append((result['log_likelihood'], tree_id))
            trees[tree_id] = result['tree']

    ml_tree_id = sorted(tree_log_likelihoods)[-1][1]
    
    with gzip.open(ml_result_file, 'w') as f:
        pickle.dump(trees[ml_tree_id], f)
        
    search_data = pd.DataFrame(tree_log_likelihoods, columns=['log_likelihood', 'tree_id'])
    search_data.to_csv(search_file, compression='gzip', sep='\t', index=False)


def annotate_posteriors(log_likelihoods_file, tree_file, nodes_file):
    ll = pd.read_csv(log_likelihoods_file, compression='gzip', sep='\t')
    tree = pickle.load(gzip.open(tree_file, 'r'))

    nodes_table = dollo.run.annotate_posteriors(ll, tree)

    nodes_table.to_csv(nodes_file, compression='gzip', sep='\t', index=False)

