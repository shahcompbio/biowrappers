from __future__ import division

from pipelines.io import make_directory

import os

from biowrappers.core import PypelinerFile, PypelinerObject

import tasks

low_mem = {'mem' : 1}
    
def single_cell_data_simulation_pipeline(sch,
                                         counts_config,
                                         tree_config,
                                         out_dir,
                                         missing_prob=0,
                                         num_counts_files=1,
                                         num_tree_files=1):
        
    dirs = {'counts' : os.path.join(out_dir, 'counts'),
            'genotypes' : os.path.join(out_dir, 'genotypes'),
            'labels' : os.path.join(out_dir, 'labels'),
            'missing' : os.path.join(out_dir, 'missing_counts'),
            'params' : os.path.join(out_dir, 'params'),
            'trees' : os.path.join(out_dir, 'trees')}
    
    for value in dirs.values():
        make_directory(value)
        
    tree_config_obj = PypelinerObject('tree_config', axes=('trees',))
    
    sch.setobj(tree_config_obj.output, get_config(tree_config, num_tree_files))
    
    base_name = get_base_name(tree_config, 'trees', 'clone-seed')

    genotypes_file = PypelinerFile(os.path.join(dirs['genotypes'], base_name + '.tsv.gz'),
                                   axes=('trees',))
    
    tree_file = PypelinerFile(os.path.join(dirs['trees'], base_name + '.nwk'),
                              axes=('trees',))
    
    sch.transform('simulate_genotypes', 
                  ('trees',), 
                  low_mem,
                  tasks.simulate_genotypes,
                  None,
                  tree_config_obj.input,
                  genotypes_file.output,
                  tree_file.output)
    
    counts_config_obj = PypelinerObject('counts_config', axes=('trees', 'counts'))

    sch.setobj(counts_config_obj.output, get_config(counts_config, num_counts_files), ('trees',))
    
    base_name = get_base_name(counts_config, 'counts', 'data-seed', prefix=base_name)

    counts_file = PypelinerFile(os.path.join(dirs['counts'], base_name + '.tsv.gz'),
                                axes=('trees', 'counts'))
    
    labels_file = PypelinerFile(os.path.join(dirs['labels'], base_name + '.tsv.gz'),
                                axes=('trees', 'counts'))
    
    params_file = PypelinerFile(os.path.join(dirs['params'], base_name + '.yaml'),
                                axes=('trees', 'counts'))
    
    sch.transform('simulate_counts', 
                  ('trees', 'counts'), 
                  low_mem,
                  tasks.simulate_count_data,
                  None,
                  counts_config_obj.input,
                  genotypes_file.input,
                  counts_file.output,
                  labels_file.output,
                  params_file.output)
    
    missing_counts_file = PypelinerFile(os.path.join(dirs['missing'], base_name + '.mp_{0}'.format(int(missing_prob * 100)) + '.tsv.gz'),
                                        axes=('trees', 'counts'))
    
    sch.transform('simulate_missing_data',
                  ('trees', 'counts'),
                  low_mem,
                  tasks.simulate_missing_data,
                  None,
                  counts_file.input,
                  missing_counts_file.output,
                  probability_missing=missing_prob,
                  seed=0)
    
    return {'counts' : counts_file,
            'genotypes' : genotypes_file,
            'labels' : labels_file,
            'missing_counts' : missing_counts_file,
            'params' : params_file,
            'tree' : tree_file}

def get_config(default_config, num_runs):
    config = {}
    
    for run_id in range(num_runs):
        config[run_id] = default_config.copy()
        
        config[run_id]['seed'] = run_id
        
    return config

def get_base_name(config, axis, axis_name, prefix=None):
    base_name = [] 
    
    for key, value in config.items():
        if type(value) is float:
            value = int(value * 100)
        
        elif type(value) is bool:
            value = int(value)
            
            continue
            
        elif type(value) is str:
            continue
        
        key = ''.join([x[0] for x in key.split('_')])
        
        base_name.append('{0}_{1}'.format(key, value))
        
    base_name.append('{0}_{1}'.format(axis_name, '{' + axis + '}'))
    
    if prefix is not None:
        base_name.insert(0, prefix)
        
    return '.'.join(base_name)