import os
import shutil
import tarfile
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

    thousand_genomes_snps = pd.read_csv(
        thousand_genomes_snps_filename, sep='\t',
        header=None, names=['chromosome', 'position', 'ref', 'alt'],
        converters={'chromosome': str})

    normal_alleles_template = os.path.join(temp_directory, normal_id + '_alleleFrequencies_chr{}.txt')
    normal_alleles_filenames = prepare_battenberg_allele_counts(
        normal_filename,
        thousand_genomes_snps,
        normal_alleles_template,
        chromosomes,
        chromosome_ids)

    tumour_alleles_template = os.path.join(temp_directory, tumour_id + '_alleleFrequencies_chr{}.txt')
    tumour_alleles_filenames = prepare_battenberg_allele_counts(
        tumour_filename,
        thousand_genomes_snps,
        tumour_alleles_template,
        chromosomes,
        chromosome_ids)

    with tarfile.open(allele_counts_filename, "w:gz") as tar:
        for filename in itertools.chain(normal_alleles_filenames, tumour_alleles_filenames):
            tar.add(filename)


def run_battenberg(
    allele_counts_filename,
    normal_id,
    tumour_id,
    results_filename,
    temp_directory,
    config):
    """ Run the battenberg method.
    """
    genome_fasta_index = config['genome_fasta_index']
    impute_info_filename = config['impute_info_filename']
    thousand_genomes_loci_directory = config['thousand_genomes_loci_directory']
    ignore_contigs_filename = config['ignore_contigs_filename']
    num_threads = config['num_threads']
    prob_loci_filename = config['prob_loci_filename']
    assembly = config['assembly']
    species = config['species']

    pypeliner.commandline.execute(
        'battenberg.pl',
        '-outdir', temp_directory,
        '-reference', genome_fasta_index,
        '-allele-counts', allele_counts_filename,
        '-gender', 'XX',
        '-impute-info', impute_info_filename,
        '-thousand-genomes-loc',  thousand_genomes_loci_directory,
        '-ignore-contigs-file', ignore_contigs_filename,
        '-t', str(num_threads),
        '-prob-loci', prob_loci_filename,
        '-tumbam', tumour_id,
        '-normbam', normal_id
        '-ra', assembly,
        '-rs', species)

    cellularity_ploidy_filename = os.path.join(temp_directory, '{}_cellularity_ploidy.txt'.format(tumour_id))
    cellularity_ploidy = pd.read_csv(cellularity_ploidy_filename, sep='\t')
    tumour_content = cellularity_ploidy['cellularity'].iloc[0]

    cn_vcf_filename = os.path.join(temp_directory, '{}_cellularity_ploidy.txt'.format(tumour_id))
    cn_data = []
    vcf_reader = vcf.Reader(gzip.open(cn_vcf_filename))
    for vcf_record in vcf_reader:
        segment_info = {}

        segment_info['chromosome'] = vcf_record.CHROM
        segment_info['start'] = vcf_record.POS
        segment_info['end'] = vcf_record.INFO['END']

        for sample in vcf_record.samples:
            if sample.sample == 'TUMOUR':
                segment_info['total_1'] = sample.data.TCN
                segment_info['minor_1'] = sample.data.MCN
                segment_info['major_1'] = segment_info['total_1'] - segment_info['minor_1']
                segment_info['fraction_1'] = sample.data.FCF

                segment_info['total_2'] = sample.data.TCS
                segment_info['minor_2'] = sample.data.MCS
                segment_info['major_2'] = segment_info['total_2'] - segment_info['minor_2']
                segment_info['fraction_2'] = sample.data.FCS

        cn_data.append(segment_info)

    cn_data = pd.DataFrame(cn_data)

    with pd.HDFStore(results_filename, 'w') as store:
        store['cn'] = cn_data
        store['mix'] = pd.Series([1. - tumour_content, tumour_content])
