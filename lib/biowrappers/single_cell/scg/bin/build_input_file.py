import yaml

def main(args):
    with open(args.template_file) as fh:
        template = yaml.load(fh)
    
    config = {'data' : {}, 
              'num_clusters' : args.num_clusters,
              'kappa_prior' : args.kappa_prior}
    
    if 'alpha_prior' in template:
        config['alpha_prior'] = template['alpha_prior']
    
    for file_name, data_type in args.in_files:
        config['data'][data_type] = template['data'][data_type]
        
        config['data'][data_type]['file']= file_name
    
    with open(args.out_file, 'w') as fh:
        yaml.dump(config, fh, default_flow_style=False)
        
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--in_files', nargs=2, action='append', required=True)
    
    parser.add_argument('--kappa_prior', required=True, type=float)
    
    parser.add_argument('--num_clusters', required=True, type=int)
    
    parser.add_argument('--template_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    args = parser.parse_args()
    
    main(args)