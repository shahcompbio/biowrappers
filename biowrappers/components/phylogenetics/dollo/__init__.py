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
    grid_search=False,
    grid_size=101,
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
    
    if grid_search:
        workflow.setobj(
            obj=pypeliner.managed.TempOutputObj('grid', 'trees', 'loss_prob'), 
            value=dict(
                zip(
                    range(grid_size), 
                    np.linspace(0, 1, grid_size),
                )
            ),
            axes=('trees',)
        ) 

        workflow.transform(
            name='compute_tree_log_likelihoods',
            axes=('trees', 'loss_prob'),
            ctx=default_ctx,
            func=tasks.compute_tree_log_likelihoods,
            args=(
                pypeliner.managed.InputFile(in_file),
                pypeliner.managed.TempInputFile('tree.pickle.gz', 'trees'),
                pypeliner.managed.TempOutputFile('log_likelhood.tsv.gz', 'trees', 'loss_prob'),
            ),
            kwargs={
                'probability_of_loss' : pypeliner.managed.TempInputObj('grid', 'trees', 'loss_prob'),
            }
        )
     
    workflow.transform(
        name='compute_optimal_tree_log_likelihoods',
        axes=('trees',),
        ctx=default_ctx,
        func=tasks.compute_tree_log_likelihoods,
        args=(
            pypeliner.managed.InputFile(in_file),
            pypeliner.managed.TempInputFile('tree.pickle.gz', 'trees'),
            pypeliner.managed.TempOutputFile('optimal_log_likelhood.tsv.gz', 'trees'),
        ),
        kwargs={
            'max_probability_of_loss' : max_probability_of_loss,
            'min_probability_of_loss' : min_probability_of_loss,
        }
    )
  
    if grid_search:
        grid_search_files = pypeliner.managed.TempInputFile('log_likelhood.tsv.gz', 'trees', 'loss_prob')
      
    else:
        grid_search_files = None
          
    workflow.transform(
        name='build_results_file', 
        axes=(), 
        ctx=default_ctx, 
        func=tasks.build_results_file, 
        args=(
            pypeliner.managed.TempInputFile('optimal_log_likelhood.tsv.gz', 'trees'),
            pypeliner.managed.OutputFile(search_file),
        ), 
        kwargs={
            'grid_search_files' : grid_search_files,
        }
    )
      
    workflow.transform(
        name='find_ml_tree',
        axes=(),
        ctx=default_ctx,
        func=tasks.find_ml_tree,
        args=(
            pypeliner.managed.InputFile(search_file),
            pypeliner.managed.TempInputFile('tree.pickle.gz', 'trees'),
            pypeliner.managed.OutputFile(tree_file)
        )
    )

    workflow.transform(
        name='get_ml_loss_prob', 
        ctx={'local' : True}, 
        func=get_loss_prob, 
        ret=pypeliner.managed.TempOutputObj('ml_loss_prob'),
        args=(
            pypeliner.managed.InputFile(search_file),
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
            '--tree_file', pypeliner.managed.InputFile(tree_file),
            '--out_file', pypeliner.managed.OutputFile(nodes_file),
            '--probability_of_loss', pypeliner.managed.TempInputObj('ml_loss_prob'),
        ),
    )
    
    return workflow

def get_loss_prob(file_name):
    
    df = pd.read_csv(file_name, compression='gzip', sep='\t')
    
    return df.iloc[0]['probability_of_loss']