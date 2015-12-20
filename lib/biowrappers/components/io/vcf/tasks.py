'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
from pypeliner.workflow import Workflow

import itertools
import pandas as pd
import pypeliner
import shutil
import time
import vcf

from biowrappers.components.utils import flatten_input

from ._merge import merge_vcfs
            
def compress_vcf(in_file, out_file):
    """ Compress a VCF file using bgzip.
    
    :param in_file: Path of uncompressed VCF file.
    :param out_file: Path were compressed VCF file will be written.
    
    """
    
    pypeliner.commandline.execute('bgzip', '-c', in_file, '>', out_file)

def filter_vcf(in_file, out_file):
    """ Filter a VCF for records with no filters set.
    
    :param in_file: Path of VCF file to filter.
    
    :param out_file: Path where filtered VCF file will be written.
    
    Note that records with the filter `PASS` will not be removed.
    
    """
    
    reader = vcf.Reader(filename=in_file)
    
    with open(out_file, 'wb') as out_fh:
        writer = vcf.Writer(out_fh, reader)
        
        for record in reader:
            if (record.FILTER is None) or (len(record.FILTER) == 0):
                writer.write_record(record)
        
        writer.close()

def finalise_vcf(in_file, compressed_file):
    """ Compress a VCF using bgzip and create index.
    
    :param workflow: pypeliner Scheduler.
    
    :param job_suffix: Suffix to add to job name. Ensures job is unique in pypeliner pipeline.
    
    :param in_file: Path of file to compressed and index.
    
    :param out_file: Path where compressed file will be written. Index file will written to `out_file` + `.tbi`.
    
    """
    
    workflow = Workflow()
    
    workflow.transform(
        name='compress_vcf',
        ctx={'mem' : 2},
        func=compress_vcf,
        args=(
            pypeliner.managed.InputFile(in_file),
            pypeliner.managed.OutputFile(compressed_file)
        )
    )
    
    workflow.transform(
        name='index_vcf',
        ctx={'mem' : 2},
        func=index_vcf,
        args=(
            pypeliner.managed.InputFile(compressed_file),
            pypeliner.managed.OutputFile(compressed_file.replace('.tmp', '') + '.tbi')
        )
    )
    
    return workflow

def index_vcf(vcf_file, index_file):
    """ Create a tabix index for a VCF file
    
    :param vcf_file: Path of VCF to create index for. Should compressed by bgzip.
    :param index_file: Path of index file.
    
    This is meant to be used from pypeliner so it does some name mangling to add .tmp to the index file.
    
    """
    
    pypeliner.commandline.execute('tabix', '-f', '-p', 'vcf', vcf_file)
    
    time.sleep(1)
    
    shutil.move(vcf_file + '.tbi', index_file)
    
    time.sleep(1)

def concatenate_vcf(in_files, out_file, variant_filter='all'):
    """ Concatenate a list of VCF files into a single file
    
    :param in_files: A dictionary of where the values are paths of VCF files to be concatenated. Files will be sorted by
                     dictionary key.
                     
    :param out_file: Path where concatenated will be written.
    
    :param variant_filter: Type of variant to include in output. Options are: `all` - all variants, `snv` - only snvs,
                           `indel` - only indels.
                           
    """

    in_files = flatten_input(in_files)
    
    reader = vcf.Reader(filename=in_files[0])
    
    with open(out_file, 'w') as out_fh:
        writer = vcf.Writer(out_fh, reader)
    
        for file_name in in_files:
            reader = vcf.Reader(filename=file_name)
            
            for record in reader:
                if variant_filter == 'all':
                    writer.write_record(record)
                
                elif variant_filter == 'indel' and record.is_indel:
                    writer.write_record(record)
                
                elif variant_filter == 'snv' and record.is_snp:
                    writer.write_record(record)
                
                else:
                    continue
                
def split_vcf(in_file, out_file_callback, lines_per_file):
    """ Split a VCF file into smaller files.
    
    :param in_file: Path of VCF file to split.
    
    :param out_file_callback: Callback function which supplies file name given index of split.
    
    :param lines_per_file: Maximum number of lines to be written per file.
    
     """
    
    def line_group(line, line_idx=itertools.count()):
        return int(next(line_idx) / lines_per_file)
    
    reader = vcf.Reader(filename=in_file)
    
    for file_idx, records in itertools.groupby(reader, key=line_group):
        file_name = out_file_callback(file_idx)
        
        with open(file_name, 'w') as out_fh:
            writer = vcf.Writer(out_fh, reader)
    
            for record in records:
                writer.write_record(record)
        
            writer.close()

def convert_vcf_to_hdf5(in_file, out_file, table_name, score_callback=None):
     
    def line_group(line, line_idx=itertools.count()):
        return int(next(line_idx) / chunk_size)
    
    chunk_size = 1000
    
    reader = vcf.Reader(filename=in_file)
    
    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    for file_idx, records in itertools.groupby(reader, key=line_group):
        df = []
        
        for record in enumerate(records):
            if score_callback is not None:
                score = score_callback(record)
            
            else:
                score = record.QUAL
            
            row = {'chrom' : record.CHROM, 'coord' : record.POS, 'score' : score}
            
            df.append(row)
        
        offset = file_idx * chunk_size
        
        df = pd.DataFrame(df, index=range(offset, offset + chunk_size))
        
        hdf_store.append(table_name, df, min_itemsize={'chrom' : 8})
    
    hdf_store.close()
    
