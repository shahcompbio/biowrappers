'''
Created on Nov 2, 2015

@author: Andrew Roth
'''
import pandas as pd
import pysam
import vcf

def get_tri_nucelotide_context(ref_genome_fasta_file, vcf_file, out_file, table_name):
    vcf_reader = vcf.Reader(filename=vcf_file)
    
    fasta_reader = pysam.Fastafile(ref_genome_fasta_file)
    
    data = []
    
    for record in vcf_reader:
        chrom = record.CHROM
        
        coord = record.POS
        
        tri_nucleotide_context = fasta_reader.fetch(chrom, coord - 2, coord + 1)
        
        data.append({'chrom' : record.CHROM, 'coord' : record.POS, 'tri_nucleotide_context' : tri_nucleotide_context})
        
    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')
    
    hdf_store[table_name] = pd.DataFrame(data, columns=['chrom', 'coord', 'tri_nucleotide_context'])
    
    hdf_store.close()