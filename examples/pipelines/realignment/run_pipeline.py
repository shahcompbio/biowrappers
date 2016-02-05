import pypeliner
import yaml

from biowrappers.pipelines.realignment import realignment_pipeline

def main(args):
    with open(args.config_file) as fh:
        config = yaml.load(fh)
    
    workflow = realignment_pipeline(
        config, 
        args.in_file, 
        args.out_file, 
    )
        
    if args.log_dir is not None:
        config['pypeliner']['tmpdir'] = args.log_dir
        
    pyp = pypeliner.app.Pypeline([], config['pypeliner'])
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--config_file', required=True)
    
    parser.add_argument('--in_file', required=True)
    
    parser.add_argument('--out_file',required=True)
    
    parser.add_argument('--log_dir', default='./')
    
    args = parser.parse_args()
    
    main(args)
