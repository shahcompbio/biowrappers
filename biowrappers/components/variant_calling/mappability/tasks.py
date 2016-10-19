'''
Created on Nov 2, 2015

@author: Andrew Roth
'''
from bx.bbi.bigwig_file import BigWigFile

import pandas as pd
import vcf

def get_mappability(
    mappability_file,
    vcf_file,
    out_file,
    table_name,
    region=None,
    append_chr=True):
    
    map_reader = BigWigFile(open(mappability_file))
    
    vcf_reader = vcf.Reader(filename=vcf_file)
    
    chrom, beg, end = _parse_region(region)
    
    if region is not None:
        vcf_reader = vcf_reader.fetch(chrom, start=beg, end=end)
    
    data = []
    
    for record in vcf_reader:
        if append_chr:
            chrom = 'chr{0}'.format(record.CHROM)
            
        else:
            chrom = record.CHROM
        
        coord = record.POS
        
        beg = coord - 100
        
        beg = max(beg, 0)
        
        end = coord + 100
        
        result = map_reader.query(chrom, beg, end, 1)
        
        if result is None:
            mappability = 0
            
        else:
            mappability = result[0]['mean']
        
        data.append({'chrom' : record.CHROM, 'coord' : record.POS, 'mappability' : mappability})
    
    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')
    
    hdf_store[table_name] = pd.DataFrame(data, columns=['chrom', 'coord', 'mappability'])
    
    hdf_store.close()


def _parse_region(region):
    if ':' not in region:
        return region, None, None

    chrom, coords = region.split(':')

    if '-' not in coords:
        return chrom, int(coords) - 1, None

    beg, end = coords.split('-')

    beg = int(beg) - 1
    end = int(end)

    return chrom, beg, end
