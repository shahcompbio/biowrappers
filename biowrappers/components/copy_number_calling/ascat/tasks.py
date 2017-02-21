import itertools
import numpy as np
import pandas as pd
import pypeliner
import remixt.seqdataio
import remixt.segalg
import remixt.analysis.haplotype
import shutil

from biowrappers.components.copy_number_calling.common.tasks import calculate_breakpoint_copy_number


def read_chromosome_lengths(chrom_info_filename):

    chrom_info = pd.read_csv(chrom_info_filename, sep='\t', compression='gzip', names=['chrom', 'length', 'twobit'])

    chrom_info['chrom'] = chrom_info['chrom'].str.replace('chr', '')

    return chrom_info.set_index('chrom')['length']


def create_segments(chrom_length, segment_length=1000):

    seg_start = np.arange(0, chrom_length, segment_length)
    seg_end = seg_start + segment_length

    segments = np.array([seg_start, seg_end]).T

    return segments


def create_segment_counts(seqdata_filename, chromosome_lengths, segment_length=1000):

    segment_counts = list()

    chromosomes = remixt.seqdataio.read_chromosomes(seqdata_filename)

    for chrom in chromosomes:
        chrom_reads = remixt.seqdataio.read_fragment_data(seqdata_filename, chromosome=chrom)

        chrom_reads.sort('start', inplace=True)

        chrom_segments = create_segments(chromosome_lengths[chrom], segment_length)

        seg_count = remixt.segalg.contained_counts(
            chrom_segments,
            chrom_reads[['start', 'end']].values,
        )

        segment_counts.append(seg_counts)

    return segment_counts


def calculate_allele_counts(seqdata_filename):

    allele_counts = list()
    
    chromosomes = remixt.seqdataio.read_chromosomes(seqdata_filename)

    for chrom in chromosomes:

        chrom_allele_counts = remixt.analysis.haplotype.read_snp_counts(seqdata_filename, chrom)
        
        chrom_allele_counts['chromosome'] = chrom

        allele_counts.append(chrom_allele_counts)

    allele_counts = pd.concat(allele_counts, ignore_index=True)

    return allele_counts


def prepare_normal_data(normal_filename, normal_logr_filename, normal_baf_filename, config):
    """ Prepare normal count and allele data
    """

    chromosome_lengths = read_chromosome_lengths(config['chrom_info_filename'])

    write_segment_count_wig(normal_wig_filename, normal_filename, chromosome_lengths, segment_length=config['window_size'])

    het_positions = infer_het_positions(normal_filename)
    het_positions.to_csv(het_positions_filename, sep='\t', index=False)


def prepare_tumour_data(tumour_filename, normal_logr_filename, normal_baf_filename, config):
    """ Prepare tumour count and allele data
    """

    chromosome_lengths = read_chromosome_lengths(config['chrom_info_filename'])

    write_segment_count_wig(tumour_wig_filename, tumour_filename, chromosome_lengths, segment_length=config['window_size'])

    het_positions = pd.read_csv(het_positions_filename, sep='\t', converters={'chromosome': str})

    tumour_allele_count = calculate_allele_counts(tumour_filename).merge(het_positions)
    write_titan_format_alleles(tumour_allele_filename, tumour_allele_count)


