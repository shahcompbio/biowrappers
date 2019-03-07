import os
import itertools
import shutil
import tarfile
import vcf
import gzip
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


def prepare_battenberg_allele_counts(
    seqdata_filename,
    thousand_genomes_snps,
    thousand_genomes_alleles_template,
    alleles_template,
    chromosomes,
    chromosome_ids):
    """ Prepare a matrix of allele counts in battenberg format.
    """

    alleles_filenames = []

    for chromosome, chromosome_id in zip(chromosomes, chromosome_ids):
        allele_counts = remixt.analysis.haplotype.read_snp_counts(seqdata_filename, chromosome)
        allele_counts['chromosome'] = chromosome

        allele_counts = allele_counts.merge(thousand_genomes_snps)
        allele_counts['ref_count'] = allele_counts['ref_count'].fillna(0).astype(int)
        allele_counts['alt_count'] = allele_counts['alt_count'].fillna(0).astype(int)

        positions = pd.read_csv(
            thousand_genomes_alleles_template.format(chromosome_id),
            sep='\t', usecols=['position'])['position'].values

        ref_counts = allele_counts.set_index(['position', 'ref'])['ref_count'].unstack(fill_value=0).reindex(index=positions, fill_value=0)
        alt_counts = allele_counts.set_index(['position', 'alt'])['alt_count'].unstack(fill_value=0).reindex(index=positions, fill_value=0)

        allele_matrix = (ref_counts + alt_counts)
        allele_matrix['Good_depth'] = allele_matrix.sum(axis=1)
        allele_matrix.reset_index(inplace=True)
        allele_matrix['chromosome'] = chromosome

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

        alleles_filename = alleles_template.format(chromosome_id)

        allele_matrix.to_csv(
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
    make_directory(temp_directory)

    chromosomes = config['chromosomes']
    chromosome_ids = config['chromosome_ids']
    thousand_genomes_snps_filename = config['thousand_genomes_snps']
    thousand_genomes_alleles_template = config['thousand_genomes_alleles_template']

    thousand_genomes_snps = pd.read_csv(
        thousand_genomes_snps_filename, sep='\t',
        header=None, names=['chromosome', 'position', 'ref', 'alt'],
        converters={'chromosome': str})

    # Ensure only SNPs are used
    thousand_genomes_snps = thousand_genomes_snps.merge(pd.DataFrame({'ref': list('ACGT')}))
    thousand_genomes_snps = thousand_genomes_snps.merge(pd.DataFrame({'alt': list('ACGT')}))

    normal_alleles_template = os.path.join(temp_directory, normal_id + '_alleleFrequencies_chr{}.txt')
    normal_alleles_filenames = prepare_battenberg_allele_counts(
        normal_filename,
        thousand_genomes_snps,
        thousand_genomes_alleles_template,
        normal_alleles_template,
        chromosomes,
        chromosome_ids)

    tumour_alleles_template = os.path.join(temp_directory, tumour_id + '_alleleFrequencies_chr{}.txt')
    tumour_alleles_filenames = prepare_battenberg_allele_counts(
        tumour_filename,
        thousand_genomes_snps,
        thousand_genomes_alleles_template,
        tumour_alleles_template,
        chromosomes,
        chromosome_ids)

    with tarfile.open(allele_counts_filename, "w:gz") as tar:
        for filename in itertools.chain(normal_alleles_filenames, tumour_alleles_filenames):
            tar.add(filename, arcname=os.path.basename(filename))


def run_battenberg(
    allele_counts_filename,
    normal_id,
    tumour_id,
    results_filename,
    temp_directory,
    config,
    **kwargs):
    """ Run the battenberg method.
    """
    make_directory(temp_directory)
    make_directory(os.path.join(temp_directory, 'tmpBattenberg', normal_id))
    make_directory(os.path.join(temp_directory, 'tmpBattenberg', tumour_id))

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
        '-normbam', normal_id,
        '-ra', assembly,
        '-rs', species)

    cellularity_ploidy_filename = os.path.join(temp_directory, '{}_cellularity_ploidy.txt'.format(tumour_id))
    cellularity_ploidy = pd.read_csv(cellularity_ploidy_filename, sep='\t')
    tumour_content = cellularity_ploidy['cellularity'].iloc[0]

    cn_vcf_filename = os.path.join(temp_directory, '{}_battenberg_cn.vcf.gz'.format(tumour_id))
    cn_data = []
    vcf_reader = vcf.Reader(filename=cn_vcf_filename)
    for vcf_record in vcf_reader:
        segment_info = {}

        segment_info['chromosome'] = vcf_record.CHROM
        segment_info['start'] = vcf_record.POS
        segment_info['end'] = vcf_record.INFO['END']

        for sample in vcf_record.samples:
            if sample.sample == 'TUMOUR':
                segment_info['major_1'] = sample.data.TCN
                segment_info['minor_1'] = sample.data.MCN
                segment_info['fraction_1'] = sample.data.FCF

                segment_info['major_2'] = sample.data.TCS
                segment_info['minor_2'] = sample.data.MCS
                segment_info['fraction_2'] = sample.data.FCS

        cn_data.append(segment_info)

    cn_data = pd.DataFrame(cn_data)

    cn_data['fraction_2'] = cn_data['fraction_2'].fillna(0.0)
    cn_data.loc[cn_data['major_2'].isnull(), 'major_2'] = cn_data.loc[cn_data['major_2'].isnull(), 'major_1']
    cn_data.loc[cn_data['minor_2'].isnull(), 'minor_2'] = cn_data.loc[cn_data['minor_2'].isnull(), 'minor_1']
    cn_data['major_2'] = cn_data['major_2'].astype(int)
    cn_data['minor_2'] = cn_data['minor_2'].astype(int)
    cn_data['total_1'] = cn_data['major_1'] + cn_data['minor_1']
    cn_data['total_2'] = cn_data['major_2'] + cn_data['minor_2']

    with pd.HDFStore(results_filename, 'w') as store:
        store['cn'] = cn_data
        store['mix'] = pd.Series([1. - tumour_content, tumour_content])
        store['brk_cn'] = pd.DataFrame(columns=['prediction_id', 'cn_1', 'cn_2'])


