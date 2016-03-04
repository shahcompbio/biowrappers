import pypeliner
import yaml

from biowrappers.pipelines.setup_reference_dbs import create_setup_reference_dbs_workflow

import biowrappers.cli as cli
    
def main(args):
    with open(args.config_file) as fh:
        # Replace {ref_path_db} in config with desired path
        config_str = fh.read()
        
        config_str = config_str.format(ref_db_path=args.ref_db_path)
        
        # Load config
        config = yaml.load(config_str)
    
    workflow = create_setup_reference_dbs_workflow(config['databases'])
    
    pyp_config = cli.load_pypeliner_config(args)
    
    pyp = pypeliner.app.Pypeline([], pyp_config)
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--config_file', required=True)
    
    parser.add_argument('--ref_db_path', required=True)
    
    cli.add_pypeliner_args(parser)
    
    args = parser.parse_args()
    
    main(args)
