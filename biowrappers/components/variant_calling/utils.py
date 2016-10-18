'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
from collections import OrderedDict

import pysam


default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']

def get_bam_regions(bam_file, split_size, chromosomes=None):
    regions = {}
    
    region_index = 0
    
    chromosome_lengths = load_bam_chromosome_lengths(bam_file, chromosomes=chromosomes)
        
    if split_size is None:
        for chrom in chromosome_lengths:
            regions[region_index] = chrom
            region_index += 1
    
    else:
        step = split_size

        for chrom, length in chromosome_lengths.iteritems():
            lside_interval = range(1, length + 1, step)
            rside_interval = range(step, length + step, step)
            
            for beg, end in zip(lside_interval, rside_interval):
                region = '{chrom}:{beg}-{end}'.format(chrom=chrom, beg=max(beg, 0), end=min(end, length))
                regions[region_index] = region
                region_index += 1
    
    return regions


def get_vcf_regions(vcf_file, split_size, chromosomes=None, zero_based=True, half_open=True):
    regions = {}
    region_index = 0
    
    chromosome_lengths = load_vcf_chromosome_lengths(vcf_file, chromosomes=chromosomes)
    
    if split_size is None:
        split_size = max(chromosome_lengths.values()) + 1
    
    for chrom, length in chromosome_lengths.iteritems():
        lside_interval = range(1, length + 1, split_size)
        rside_interval = range(split_size, length + split_size, split_size)
        
        for beg, end in zip(lside_interval, rside_interval):
            end = min(end, length)
            
            if zero_based:
                beg -= 1
                end -= 1
            
            if half_open:
                end += 1
            
            regions[region_index] = (chrom, beg, end)
            region_index += 1

    return regions


def load_vcf_chromosome_lengths(file_name, chromosomes=None):
    chromosome_lengths = OrderedDict()

    vcf_reader = vcf.Reader(filename=vcf_file, compressed=True)

    if chromosomes is None:
        chromosomes = vcf_reader.contigs.keys()

    else:
        chromosomes = chromosomes

    for chrom, contig in vcf_reader.contigs.iteritems():
        assert chrom == contig.id

        if chrom not in chromosomes:
            continue

        chromosome_lengths[chrom] = int(contig.length)

    return chromosome_lengths


def load_bam_chromosome_lengths(file_name, chromosomes=None):
    chromosome_lengths = OrderedDict()
    
    bam = pysam.Samfile(file_name, 'rb')

    if chromosomes is None:
        chromosomes = bam.references
    
    else:
        chromosomes = chromosomes
    
    for chrom, length in zip(bam.references, bam.lengths):
        if chrom not in chromosomes:
            continue
        
        chromosome_lengths[chrom] = int(length)
    
    return chromosome_lengths
