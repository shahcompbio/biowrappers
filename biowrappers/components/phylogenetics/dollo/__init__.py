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
        name='create_tree_groups', 
        ctx=default_ctx, 
        func=tasks.create_tree_groups, 
        args=(
            pypeliner.managed.InputFile(in_file),
            pypeliner.managed.TempOutputFile('trees.pickle.gz', 'tree_groups'),
        ),
    )
    
    workflow.transform(
        name='compute_optimal_tree_log_likelihoods',
        axes=('tree_groups',),
        ctx=default_ctx,
        func=tasks.compute_tree_log_likelihoods,
        args=(
            pypeliner.managed.InputFile(in_file),
            pypeliner.managed.TempInputFile('trees.pickle.gz', 'tree_groups'),
            pypeliner.managed.TempOutputFile('results.pickle.gz', 'tree_groups'),
        ),
        kwargs={
            'max_probability_of_loss' : max_probability_of_loss,
            'min_probability_of_loss' : min_probability_of_loss,
        }
    )
  
    workflow.transform(
        name='select_ml_tree',
        ctx=default_ctx,
        func=tasks.select_ml_tree,
        args=(
            pypeliner.managed.TempInputFile('results.pickle.gz', 'tree_groups'),
            pypeliner.managed.OutputFile(tree_file),
            pypeliner.managed.OutputFile(search_file),
        )
    )

    workflow.transform(
        name='annotate_posteriors',
        ctx=default_ctx,
        func=tasks.annotate_posteriors,
        args=(
            pypeliner.managed.InputFile(in_file),
            pypeliner.managed.InputFile(tree_file),
            pypeliner.managed.OutputFile(nodes_file),
        ),
    )
    
    return workflow

