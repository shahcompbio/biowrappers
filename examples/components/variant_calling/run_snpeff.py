from pypeliner.workflow import Workflow

import pypeliner

from biowrappers.components.io.hdf5.tasks import convert_hdf5_to_tsv

import biowrappers.cli as cli
import biowrappers.components.variant_calling.snpeff as snpeff

def main(args):
    config = cli.load_pypeliner_config(args)
    
    pyp = pypeliner.app.Pypeline([], config)
    
    workflow = Workflow()
    
    workflow.subworkflow(
        name='snpeff',
        func=snpeff.create_snpeff_annotation_workflow, 
        args=(
              pypeliner.managed.InputFile(args.target_vcf_file),
              pypeliner.managed.TempOutputFile('snpeff.h5')
        ), 
        kwargs={
            'data_base' : args.data_base,
            'split_size' : args.split_size,
            'table_name' : 'snpeff'
        }
    )

    workflow.transform(
        name='convert_to_tsv', 
        func=convert_hdf5_to_tsv,
        ctx={'mem' : 2},
        args=(
            pypeliner.managed.TempInputFile('snpeff.h5'),
            'snpeff',
            pypeliner.managed.OutputFile(args.out_file)
        ),
        kwargs={
            'compress' : True,
            'index' : False
        }                               
    )
    
    pyp.run(workflow)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--target_vcf_file', required=True)
    
    parser.add_argument('--out_file', required=True)
    
    parser.add_argument('--data_base', default='GRCh37.75')
    
    parser.add_argument('--split_size', default=int(1e3), type=int)
    
    cli.add_pypeliner_args(parser)
    
    args = parser.parse_args()
    
    main(args)
    
