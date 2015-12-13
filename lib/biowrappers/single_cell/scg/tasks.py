'''
Created on Oct 30, 2015

@author: Andrew Roth
'''
import gzip
import pandas as pd
import yaml

def write_scg_config_file(in_file, out_file, template_file, num_clusters=40, kappa_prior=1):
    with open(template_file) as fh:
        template = yaml.load(fh)
    
    config = {'data' : {},
              'num_clusters' : num_clusters,
              'kappa_prior' : kappa_prior}
    
    if 'alpha_prior' in template:
        config['alpha_prior'] = template['alpha_prior']
    
#     for data_type, file_name in zip(data_types, in_files):
#         config['data'][data_type] = template['data'][data_type]
#         
#         config['data'][data_type]['file'] = file_name
#         
    config['data']['snv'] = template['data']['snv']
    
    config['data']['snv']['file'] = in_file
    
    with open(out_file, 'w') as fh:
        yaml.dump(config, fh, default_flow_style=False)
        
def write_scg_samples_file(in_files, out_file):
    cells = []
    
    for file_name in in_files:
        df = pd.read_csv(file_name, compression='gzip', sep='\t')
        
        cells.append(set(df['cell_id']))
    
    cells = list(set.intersection(*cells))

    data = []
    
    for cell_id in cells:
        data.append({'cell_id' : cell_id, 'sample' : cell_id.split(':')[0]})
    
    data = pd.DataFrame(data)
    
    with gzip.GzipFile(out_file, 'w') as fh:
        data.to_csv(fh, index=False, sep='\t')
