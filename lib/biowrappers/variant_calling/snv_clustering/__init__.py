import os
import pypeliner
import yaml

from biowrappers.variant_calling.snv_allele_counts import snv_allele_counts_pipeline
from biowrappers.variant_calling.utils import default_chromosomes

import tasks
from numpy.random.mtrand import seed

def snv_allele_count_clustering_pipeline(
    sch,
    counts_file,
    out_file,
    num_restarts=10):
    
    sch.transform(
        'build_input',
        (),
        {'mem' : 4},
        tasks.build_input_file,
        None,
        pypeliner.managed.InputFile(counts_file),
        pypeliner.managed.TempOutputFile('input.h5')
    )
  
    sch.setobj(
        pypeliner.managed.TempOutputObj('seeds_obj', 'seeds'),
        get_restart_seeds(num_restarts)
    )
    
    sch.transform(
        'seed_search',
        ('seeds',),
        {'mem' : 4},
        tasks.cluster_snv_counts,
        None,
        pypeliner.managed.TempInputFile('input.h5'),
        lower_bound_file=pypeliner.managed.TempOutputFile('lower_bound.yaml', 'seeds'),
        seed=pypeliner.managed.TempInputObj('seeds_obj', 'seeds')
    )
    
    sch.transform(
        'find_best_seed',
        (),
        {'mem' : 2},
        find_best_run,
        pypeliner.managed.TempOutputObj('best_seed'),
        pypeliner.managed.TempInputFile('lower_bound.yaml', 'seeds')
    )
    
    sch.transform(
        'write_best_run',
        (),
        {'mem' : 4},
        tasks.cluster_snv_counts,
        None,
        pypeliner.managed.TempInputFile('input.h5'),
        params_file=pypeliner.managed.OutputFile(out_file),
        seed=pypeliner.managed.TempInputObj('best_seed')
    )

def get_restart_seeds(num_restarts):
    
    seeds = {}
    
    for i in range(num_restarts):
        seeds[i] = i
    
    return seeds

def find_best_run(lower_bounds_files):
    
    best_lower_bound = float('-inf')
    
    best_seed = None
    
    for seed, file_name in lower_bounds_files.items():
        with open(file_name) as fh:
            params = yaml.load(fh)
            
            if params['lower_bound'] > best_lower_bound:
                best_lower_bound = params['lower_bound']
                
                best_seed = seed
    
    return best_seed