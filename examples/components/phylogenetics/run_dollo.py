import pypeliner

import biowrappers.cli as cli
import biowrappers.components.phylogenetics.dollo as dollo

def main(args):
    config = cli.load_pypeliner_config(args)
    
    pyp = pypeliner.app.Pypeline([], config)
    
    workflow = dollo.get_tree_search_workflow(
        args.in_file, 
        args.search_file,
        args.tree_file, 
        grid_search=args.grid_search, 
        grid_size=args.grid_size, 
        max_probability_of_loss=args.max_probability_of_loss, 
        min_probability_of_loss=args.min_probability_of_loss
    )
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--in_file', required=True)
    
    parser.add_argument('--search_file', required=True)
    
    parser.add_argument('--tree_file', required=True)
    
    parser.add_argument('--grid_search', default=False, action='store_true')
    
    parser.add_argument('--grid_size', default=101, type=int)
    
    parser.add_argument('--max_probability_of_loss', type=float, default=0.5)

    parser.add_argument('--min_probability_of_loss', type=float, default=0.0)
    
    parser.add_argument('--num_probability_of_loss', type=int, default=51)

    cli.add_pypeliner_args(parser)
    
    args = parser.parse_args()
    
    main(args)
