'''
Created on Nov 2, 2015

@author: Andrew Roth
'''

import pandas as pd
import pypeliner
import vcf

def run_snpeff(db, in_vcf_file, out_file):
    
    cmd = [
        'snpEff',
        '-noStats',
        '-noLog',
        '-classic',
        '-Xms1g',
        '-Xmx12g',
        db,
        in_vcf_file,
        '>',
        out_file
    ]
    
    pypeliner.commandline.execute(*cmd)
    
#=======================================================================================================================
# TSV output
#=======================================================================================================================
effect_keys = (
               'effect_impact',
               'functional_class',
               'codon_change',
               'amino_acid_change',
               'amino_acid_length',
               'gene_name',
               'biotype',
               'gene_coding',
               'feature_id',
               'exon_rank',
               'genotype_number'
               )
                          
effect_impact = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')

cols = ['chrom', 'coord', 'ref', 'alt', 'effect'] + list(effect_keys)

def convert_vcf_to_table(in_file, out_file, table_name):
    data = []
    
    for out_row in _get_annotations_table(in_file):
        data.append(out_row)
    
    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')
    
    hdf_store[table_name] = pd.DataFrame(data, columns=cols)
    
    hdf_store.close()
    
def _get_annotations_table(file_name):
    reader = vcf.VCFReader(filename=file_name)

    for record in reader:
        if 'EFF' not in record.INFO:
            continue
        
        for status, values in _parse_snpeff_field(record.INFO['EFF']):
            out_row = {
                       'chrom' : record.CHROM,
                       'coord' : record.POS,
                       'ref' : record.REF,
                       'alt' : ','.join([str(x) for x in record.ALT]),
                       'effect' : status,
                       }
            
            for key in effect_keys:
                out_row[key] = values[key]
        
            yield out_row

def _parse_snpeff_field(field):
    for effect_field in field:
        effect_status, values = effect_field.split('(')
        
        values = values.replace(')', '')
        
        values = values.split('|')
        
        effect_values = dict(zip(effect_keys, values))

        yield effect_status, effect_values
