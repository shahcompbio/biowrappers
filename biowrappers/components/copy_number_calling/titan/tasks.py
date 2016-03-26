import itertools
import numpy as np
import pandas as pd
import pypeliner
import remixt.seqdataio
import remixt.segalg
import remixt.analysis.haplotype
import shutil

def read_chromosome_lengths(chrom_info_filename):

    chrom_info = pd.read_csv(chrom_info_filename, sep='\t', compression='gzip', names=['chrom', 'length', 'twobit'])

    chrom_info['chrom'] = chrom_info['chrom'].str.replace('chr', '')

    return chrom_info.set_index('chrom')['length']


def create_segments(chrom_length, segment_length=1000):

    seg_start = np.arange(0, chrom_length, segment_length)
    seg_end = seg_start + segment_length

    segments = np.array([seg_start, seg_end]).T

    return segments


def write_segment_count_wig(wig_filename, seqdata_filename, chromosome_lengths, segment_length=1000):

    with open(wig_filename, 'w') as wig:

        chromosomes = remixt.seqdataio.read_chromosomes(seqdata_filename)

        for chrom in chromosomes:

            wig.write('fixedStep chrom={0} start=1 step={1} span={1}\n'.format(chrom, segment_length))

            chrom_reads = remixt.seqdataio.read_filtered_fragment_data(seqdata_filename, chromosome=chrom)

            chrom_reads.sort('start', inplace=True)

            chrom_segments = create_segments(chromosome_lengths[chrom], segment_length)

            seg_count = remixt.segalg.contained_counts(
                chrom_segments,
                chrom_reads[['start', 'end']].values,
            )

            wig.write('\n'.join([str(c) for c in seg_count]))
            wig.write('\n')


def calculate_allele_counts(seqdata_filename):

    allele_counts = list()
    
    chromosomes = remixt.seqdataio.read_chromosomes(seqdata_filename)

    for chrom in chromosomes:

        chrom_allele_counts = remixt.analysis.haplotype.read_snp_counts(seqdata_filename, chrom)
        
        chrom_allele_counts['chromosome'] = chrom

        allele_counts.append(chrom_allele_counts)

    allele_counts = pd.concat(allele_counts, ignore_index=True)

    return allele_counts


def infer_het_positions(seqdata_filename):

    allele_count = calculate_allele_counts(seqdata_filename)
    
    remixt.analysis.haplotype.infer_snp_genotype(allele_count)

    het_positions = allele_count.loc[allele_count['AB'] == 1, ['chromosome', 'position']]

    return het_positions


def write_titan_format_alleles(allele_filename, allele_count):

    allele_count = allele_count.rename(columns={
        'chromosome': 'chr',
        'ref_count': 'refCount',
        'alt_count': 'NrefCount',
    })

    allele_count['refBase'] = 'A'
    allele_count['NrefBase'] = 'T'

    allele_count = allele_count[[
        'chr',
        'position',
        'refBase',
        'refCount',
        'NrefBase',
        'NrefCount',
    ]]

    allele_count['refCount'] = allele_count['refCount'].astype(int)
    allele_count['NrefCount'] = allele_count['NrefCount'].astype(int)

    allele_count.to_csv(allele_filename, sep='\t', index=False, header=False)


def read_titan_params(params_filename):

    params = dict()

    with open(params_filename, 'r') as params_file:
        for line in params_file:
            key, value = line.split(':')
            params[key] = np.array(value.split()).astype(float)

    return params


def prepare_normal_data(normal_filename, normal_wig_filename, het_positions_filename, config):
    """ Prepare normal count data and infer het positions
    """

    chromosome_lengths = read_chromosome_lengths(config['chrom_info_filename'])

    write_segment_count_wig(normal_wig_filename, normal_filename, chromosome_lengths, segment_length=config['window_size'])

    het_positions = infer_het_positions(normal_filename)
    het_positions.to_csv(het_positions_filename, sep='\t', index=False)


def prepare_tumour_data(tumour_filename, het_positions_filename, tumour_wig_filename, tumour_allele_filename, config):
    """ Prepare tumour count and allele data
    """

    chromosome_lengths = read_chromosome_lengths(config['chrom_info_filename'])

    write_segment_count_wig(tumour_wig_filename, tumour_filename, chromosome_lengths, segment_length=config['window_size'])

    het_positions = pd.read_csv(het_positions_filename, sep='\t')

    tumour_allele_count = calculate_allele_counts(tumour_filename).merge(het_positions)
    write_titan_format_alleles(tumour_allele_filename, tumour_allele_count)


def create_intialization_parameters(config):
    """ Initialize parameter sweep
    """

    normal_contamination = config.get('normal_contamination', [0.2, 0.4, 0.6, 0.8])
    
    num_clusters = config.get('num_clusters', [1, 2, 3, 4, 5])
    
    ploidy = config.get('ploidy', [1, 2, 3, 4])

    init_param_values = itertools.product(
        normal_contamination,
        num_clusters,
        ploidy,
    )

    init_param_cols = [
        'normal_contamination',
        'num_clusters',
        'ploidy',
    ]

    init_params = {}
    for idx, params in enumerate(init_param_values):
        init_params[idx] = dict(zip(init_param_cols, params))

    return init_params


def run_titan(init_params, normal_wig_filename, tumour_wig_filename, tumour_allele_filename, cn_filename, params_filename, config):
    """ Run the analysis with specific initialization parameters

    """

    titan_cmd = [
        'run_titan.R',
        tumour_allele_filename,
        normal_wig_filename,
        tumour_wig_filename,
        cn_filename,
        params_filename,
        '--estimate_clonal_prevalence',
        '--estimate_normal_contamination',
        '--estimate_ploidy',
        '--max_copy_number', str(config.get('max_copy_number', 5)),
        '--max_depth', str(config.get('max_depth', int(1e4))),
        '--normal_contamination', str(init_params['normal_contamination']),
        '--num_clusters', str(init_params['num_clusters']),
        '--ploidy', str(init_params['ploidy']),
        '--gc_wig_file', config['gc_wig'],
        '--map_wig_file', config['mappability_wig'],
    ]

    pypeliner.commandline.execute(*titan_cmd)


def select_solution(
    init_params,
    cn_filename,
    params_filename,
    results_filename,
    titan_loci_filename,
    titan_segments_filename,
    titan_igv_filename,
    config,
    sample_id,
):
    """ Select optimal copy number and mixture

    """

    init_params = pd.DataFrame.from_dict(init_params, orient='index')

    for init_param_idx, row in init_params.iterrows():
        titan_params_filename = params_filename[init_param_idx]

        titan_params = read_titan_params(titan_params_filename)

        init_params.loc[init_param_idx, 'model_selection_index'] = titan_params['S_Dbw validity index (Both)'][0]
        init_params.loc[init_param_idx, 'norm_contam_est'] = titan_params['Normal contamination estimate'][0]

        cell_prev_field = 'Clonal cluster cellular prevalence Z={0}'.format(int(row['num_clusters']))
        for idx, cell_prev in enumerate(titan_params[cell_prev_field]):
            init_params.loc[init_param_idx, 'cell_prev_est_{0}'.format(idx + 1)] = cell_prev

    best_idx = init_params['model_selection_index'].argmin()

    mix = [init_params.loc[best_idx, 'norm_contam_est'], ]
    
    for i in range(int(init_params.loc[best_idx, 'num_clusters'])):
        mix.append(init_params.loc[best_idx, 'cell_prev_est_{0}'.format(i + 1)])
    
    shutil.copyfile(cn_filename[best_idx], titan_loci_filename)
    
    pypeliner.commandline.execute(
        'createTITANsegmentfiles.pl',
        '-id', sample_id,
        '-i', cn_filename[best_idx],
        '-o', titan_segments_filename,
        '-igv', titan_igv_filename,
    )

    cn_data = pd.read_csv(
        titan_segments_filename,
        sep='\t', converters={'Chromosome': str}
    )

    cn_columns = {
        'Chromosome': 'chromosome',
        'Start_Position(bp)': 'start',
        'End_Position(bp)': 'end',
        'Copy_Number': 'total_1',
        'MajorCN': 'major_1',
        'MinorCN': 'minor_1',
        'Clonal_Cluster': 'clone',
        'Clonal_Frequency' : 'prevalence'
    }

    cn_data = cn_data.rename(columns=cn_columns)[cn_columns.values()]

    cn_data['clone'] = cn_data['clone'].fillna(-1).astype(int)

    cn_data['prevalence'] = cn_data['prevalence'].fillna(0).astype(float)

    with pd.HDFStore(results_filename, 'w') as store:
        store['mix'] = pd.Series(mix)
        store['cn'] = cn_data

