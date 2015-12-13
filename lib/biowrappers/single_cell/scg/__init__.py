from pipelines.io import get_ancestor_directory

import os

from biowrappers.core import PypelinerFile

import tasks

low_mem = {'mem' : 1}

med_mem = {'mem' : 4}

def single_cell_genotyper_analysis_pipeline(sch, in_file, template_file):
    
    out_dir = get_ancestor_directory(in_file.path, level=2)
    
    config_file = os.path.join(out_dir, 'analysis', 'config', os.path.basename(in_file.path))
    
    config_file = config_file.replace('.tsv.gz', '.yaml')
    
    config_file = PypelinerFile(config_file, in_file.axes)
    
    sch.transform('write_config_file', 
                  in_file.axes, 
                  low_mem,
                  tasks.write_scg_config_file,
                  None,
                  in_file.input,
                  config_file.output,
                  template_file)