import os
import shutil
import numpy as np
import pandas as pd

import pypeliner

import remixt.seqdataio
import remixt.segalg
import remixt.analysis.haplotype
import remixt.analysis.experiment
import remixt.cn_model

from biowrappers.components.utils import make_directory
from biowrappers.components.copy_number_calling.common.tasks import calculate_breakpoint_copy_number
from biowrappers.components.copy_number_calling.common.utils import calculate_allele_counts


thousand_genomes_loci = '{ref_data_dir}/data/ref_db/battenberg/1000genomesloci/1000genomesloci2012_chr{chromosome}.txt'
thousand_genomes_snps =  '{ref_data_dir}/data/ref_db/remixt/thousand_genomes_snps.tsv'


def prepare_battenberg_allele_counts(
    seqdata_filename,
    thousand_genomes_snps,
    alleles_template,
    chromosomes,
    chromosome_ids):
    """ Prepare a matrix of allele counts in battenberg format.
    """

    allele_counts = calculate_allele_counts(seqdata_filename)
    allele_counts = normal_counts.merge(thousand_genomes_snps)

    ref_counts = allele_counts.set_index(['chromosome', 'position', 'ref'])['ref_counts'].unstack(fill_value=0)
    alt_counts = allele_counts.set_index(['chromosome', 'position', 'alt'])['alt_counts'].unstack(fill_value=0)

    allele_matrix = (ref_counts + alt_counts)
    allele_matrix['Good_depth'] = allele_matrix.sum(axis=1)
    allele_matrix.reset_index(inplace=True)
    allele_matrix.rename(
        columns={
            'chromosome': '#CHR',
            'position': 'POS',
            'A': 'Count_A',
            'C': 'Count_C',
            'G': 'Count_G',
            'T': 'Count_T',
        },
        inplace=True,
    )

    alleles_filenames = []

    for chromosome, chromosome_id in zip(chromosomes, chromosome_ids):
        alleles_filename = alleles_template.format(chromosome_id)
        chr_allele_matrix = allele_matrix.loc[allele_matrix['#CHR'] == chromosome]
        chr_allele_matrix.to_csv(
            alleles_filename,
            sep='\t', index=False,
            columns=['#CHR', 'POS', 'Count_A', 'Count_C', 'Count_G', 'Count_T', 'Good_depth'])
        alleles_filenames.append(alleles_filename)

    return alleles_filenames


def prepare_data(
    normal_filename,
    tumour_filename,
    normal_id,
    tumour_id,
    allele_counts_filename,
    temp_directory,
    config):
    """ Create battenberg input data.
    """

    chromosomes = config['chromosomes']
    chromosome_ids = config['chromosome_ids']
    thousand_genomes_snps_filename = config['thousand_genomes_snps']
    alleles_template = os.path.join(temp_directory, config['alleles_template'])

    thousand_genomes_snps = pd.read_csv(
        thousand_genomes_snps_filename, sep='\t',
        header=None, names=['chromosome', 'position', 'ref', 'alt'],
        converters={'chromosome': str})

    normal_alleles_filenames = prepare_battenberg_allele_counts(
        normal_filename,
        thousand_genomes_snps,
        alleles_template,
        chromosomes,
        chromosome_ids)

    tumour_alleles_filenames = prepare_battenberg_allele_counts(
        tumour_filename,
        thousand_genomes_snps,
        alleles_template,
        chromosomes,
        chromosome_ids)

    



def run_battenberg(
    allele_counts_filename,
    normal_id,
    sample_id,
    results_filename,
    temp_directory,
    config):
    """ Run the battenberg method.
    """
    
