import pypeliner
import yaml

from biowrappers.pipelines.setup_reference_dbs import create_setup_reference_dbs_workflow

import biowrappers.cli as cli
    
def main(args):
    with open(args.config_file) as fh:
        config = yaml.load(fh)
    
    workflow = create_setup_reference_dbs_workflow(config['databases'])
    
    pyp_config = cli.load_pypeliner_config(args)
    
    pyp = pypeliner.app.Pypeline([], pyp_config)
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--config_file', required=True)
    
    cli.add_pypeliner_args(parser)
    
    args = parser.parse_args()
    
    main(args)
