import os
import numpy as np
import pandas as pd

import pypeliner.commandline
import remixt.seqdataio
import remixt.analysis.haplotype

from biowrappers.components.copy_number_calling.common.tasks import calculate_breakpoint_copy_number
import biowrappers.components.utils as utils


def calculate_allele_counts(seqdata_filename, chromosomes=None):
    allele_counts = list()

    if chromosomes is None:
        chromosomes = remixt.seqdataio.read_chromosomes(seqdata_filename)

    for chrom in chromosomes:
        chrom_allele_counts = remixt.analysis.haplotype.read_snp_counts(seqdata_filename, chrom)
        chrom_allele_counts['chromosome'] = chrom
        allele_counts.append(chrom_allele_counts)

    allele_counts = pd.concat(allele_counts, ignore_index=True)

    return allele_counts


def write_theta_format_alleles(allele_filename, allele_count):
    allele_count = allele_count[[
        'chromosome',
        'position',
        'ref_count',
        'alt_count',
    ]]
    
    allele_count['ref_count'] = allele_count['ref_count'].astype(int)
    allele_count['alt_count'] = allele_count['alt_count'].astype(int)

    allele_count.to_csv(allele_filename, sep='\t', index=False, header=False)


bicseq2_reads = 'chr{}.seq'
bicseq2_normed = 'chr{}.norm.bin'


def run_bicseq2_norm(prefix, seqdata_filename, config, tmp_directory):
    read_length = config['read_length']

    fragment_means = []
    counts = []
    for chromosome in config['chromosomes']:
        chrom_reads = remixt.seqdataio.read_filtered_fragment_data(seqdata_filename, chromosome=chromosome)

        fragment_means.append((chrom_reads['end'] - chrom_reads['start']).mean())
        counts.append(len(chrom_reads.index))

        chrom_reads[['start']].to_csv(prefix + bicseq2_reads.format(chromosome), index=False, header=False)

    fragment_means = np.array(fragment_means)
    counts = np.array(counts)

    fragment_length = int(np.sum(fragment_means * counts) / counts.sum())

    norm_config_filename = os.path.join(tmp_directory, 'norm.config')

    with open(norm_config_filename, 'w') as f:
        f.write('chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n')
        for chromosome in config['chromosomes']:
            f.write('{}\t'.format(chromosome))
            f.write('{}\t'.format(config['chromosome_template'].format(chromosome)))
            f.write('{}\t'.format(config['mappability_template'].format(chromosome)))
            f.write('{}\t'.format(prefix + bicseq2_reads.format(chromosome)))
            f.write('{}\n'.format(prefix + bicseq2_normed.format(chromosome)))

    norm_output_filename = prefix + 'norm.output'

    pypeliner.commandline.execute(
        'bicseq2-norm',
        '-l', read_length,
        '-s', fragment_length,
        '-b', str(config['bin_size']),
        '--tmp', tmp_directory,
        norm_config_filename,
        norm_output_filename,
    )


def run_bicseq2_seg(seg_output_filename, normal_filename, tumour_filename, config, tmp_directory):
    utils.make_directory(tmp_directory)

    normal_prefix = os.path.join(tmp_directory, 'normal.')
    tumour_prefix = os.path.join(tmp_directory, 'tumour.')

    run_bicseq2_norm(normal_prefix, normal_filename, config, tmp_directory)
    run_bicseq2_norm(tumour_prefix, tumour_filename, config, tmp_directory)

    seg_config_filename = os.path.join(tmp_directory, 'seg.config')

    with open(seg_config_filename, 'w') as f:
        f.write('chromName\tbinFileNorm.Case\tbinFileNorm.Control\n')
        for chromosome in config['chromosomes']:
            f.write('{}\t'.format(chromosome))
            f.write('{}\t'.format(tumour_prefix + bicseq2_normed.format(chromosome)))
            f.write('{}\n'.format(normal_prefix + bicseq2_normed.format(chromosome)))

    pypeliner.commandline.execute(
        'bicseq2-seg',
        '--control',
        '--tmp', tmp_directory,
        seg_config_filename,
        seg_output_filename,
    )


def write_results(theta_prefix, output_filename, **kwargs):
    breakpoints_filename = kwargs.get('breakpoints_filename', None)
    num_clones = kwargs.get('num_clones', None)

    store = pd.HDFStore(output_filename, 'w')

    solution_name = 'BEST'
    if num_clones is not None:
        solution_name = 'n{}'.format(num_clones)

    theta2_results_filename = '.'.join([theta_prefix, solution_name, 'results'])
    theta2_results = pd.read_csv(theta2_results_filename, sep='\t').rename(columns={'#NLL':'NLL'})

    store['full'] = theta2_results

    best_idx = theta2_results['NLL'].argmin()

    best_frac = theta2_results.loc[best_idx, 'mu']
    best_frac = best_frac.split(',')

    store['mix'] = pd.Series(np.array(best_frac))

    best_cn = theta2_results.loc[best_idx, 'C']
    best_cn = [a.split(',') for a in best_cn.split(':')]
    best_cn = np.array(best_cn).astype(int).T

    theta2_seg_filename = theta_prefix + '.n2.withBounds'
    cn_data = pd.read_csv(theta2_seg_filename, sep='\t')
    cn_data.rename(columns={'chrm': 'chromosome'}, inplace=True)
    cn_data['chromosome'] = cn_data['chromosome'].astype(str)

    for m in range(best_cn.shape[0]):
        cn_data['total_{}'.format(m + 1)] = best_cn[m]

    store['cn'] = cn_data

    if breakpoints_filename is not None:
        store['brk_cn'] = calculate_breakpoint_copy_number(breakpoints_filename, cn_data)


def run_theta(output_filename, normal_filename, tumour_filename, bicseq2_seg_filename, config, tmp_directory, **kwargs):
    utils.make_directory(tmp_directory)

    normal_allele_filename = os.path.join(tmp_directory, 'normal_alleles.tsv')
    tumour_allele_filename = os.path.join(tmp_directory, 'tumour_alleles.tsv')

    normal_allele_count = calculate_allele_counts(normal_filename, chromosomes=config['chromosomes'])
    tumour_allele_count = calculate_allele_counts(tumour_filename, chromosomes=config['chromosomes'])

    positions = pd.merge(
        normal_allele_count[['chromosome', 'position']].drop_duplicates(),
        tumour_allele_count[['chromosome', 'position']].drop_duplicates(),
        how='inner')

    normal_allele_count = normal_allele_count.merge(positions, how='inner')
    tumour_allele_count = tumour_allele_count.merge(positions, how='inner')

    write_theta_format_alleles(normal_allele_filename, normal_allele_count)
    write_theta_format_alleles(tumour_allele_filename, tumour_allele_count)

    theta_seg_filename = os.path.join(tmp_directory, 'theta_input.seg')
    pypeliner.commandline.execute(
        'BICSeqToTHetA',
        bicseq2_seg_filename,
        '-OUTPUT_PREFIX', theta_seg_filename,
    )

    theta_prefix = os.path.join(tmp_directory, 'theta_results')
    pypeliner.commandline.execute(
        'RunTHetA', '--FORCE',
        os.path.abspath(theta_seg_filename + '.all_processed'),
        '--TUMOR_FILE', os.path.abspath(normal_allele_filename),
        '--NORMAL_FILE', os.path.abspath(tumour_allele_filename),
        '--DIR', os.path.abspath(tmp_directory),
        '--OUTPUT_PREFIX', 'theta_results',
    )

    write_results(theta_prefix, output_filename, **kwargs)

