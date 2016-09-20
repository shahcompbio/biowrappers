from pypeliner.workflow import Workflow

import numpy as np
import pandas as pd
import pypeliner

import tasks

default_ctx = {'mem' : 2, 'num_retry' : 3, 'mem_retry_increment' : 2}

def create_tree_search_workflow(
    in_file,
    nodes_file,
    search_file,
    tree_file,
    max_probability_of_loss=0.5,
    min_probability_of_loss=0.0):

    workflow = Workflow()
    
    workflow.transform(
        name='create_trees', 
        axes=(), 
        ctx=default_ctx, 
        func=tasks.create_trees, 
        args=(
            pypeliner.managed.InputFile(in_file),
            pypeliner.managed.TempOutputFile('tree.pickle.gz', 'trees'),
            pypeliner.managed.TempSpace('trees_tmp'),
        ),
    )
    
    workflow.transform(
        name='compute_optimal_tree_log_likelihoods',
        axes=('trees',),
        ctx=default_ctx,
        func=tasks.compute_tree_log_likelihoods,
        args=(
            pypeliner.managed.InputFile(in_file),
            pypeliner.managed.TempInputFile('tree.pickle.gz', 'trees'),
            pypeliner.managed.TempOutputFile('solution.pickle.gz', 'trees'),
        ),
        kwargs={
            'max_probability_of_loss' : max_probability_of_loss,
            'min_probability_of_loss' : min_probability_of_loss,
        }
    )
  
    workflow.transform(
        name='select_ml_tree',
        axes=(),
        ctx=default_ctx,
        func=tasks.select_ml_tree,
        args=(
            pypeliner.managed.TempInputFile('solution.pickle.gz', 'trees'),
            pypeliner.managed.TempOutputFile('ml_solution.pickle.gz'),
        )
    )

    workflow.commandline(
        name='annotate_posterirors',
        axes=(),
        ctx=default_ctx,
        args=(
            'PyDollo',
            'annotate_posteriors',
            '--log_likelihoods_file', pypeliner.managed.InputFile(in_file),
            '--solution_file', pypeliner.managed.TempInputFile('ml_solution.pickle.gz'),
            '--out_file', pypeliner.managed.OutputFile(nodes_file),
        ),
    )
    
    return workflow

def get_loss_prob(file_name):
    
    df = pd.read_csv(file_name, compression='gzip', sep='\t')
    
    return df.iloc[0]['probability_of_loss']