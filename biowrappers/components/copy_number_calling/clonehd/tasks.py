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


def read_chromosome_lengths(chrom_info_filename):

    chrom_info = pd.read_csv(chrom_info_filename, sep='\t', compression='gzip', names=['chrom', 'length', 'twobit'])

    chrom_info['chrom'] = chrom_info['chrom'].str.replace('chr', '')

    return chrom_info.set_index('chrom')['length']


def create_segments(chrom_length, segment_length=1000):

    seg_start = np.arange(0, chrom_length, segment_length)
    seg_end = seg_start + segment_length

    segments = pd.DataFrame({'start': seg_start, 'end': seg_end})

    return segments


def write_cna(cna_filename, seqdata_filename, chromosome_lengths, segment_length=1000):

    with open(cna_filename, 'w') as cna:

        chromosomes = remixt.seqdataio.read_chromosomes(seqdata_filename)

        for chrom in chromosomes:

            reads = remixt.seqdataio.read_fragment_data(seqdata_filename, chromosome=chrom)

            reads.sort('start', inplace=True)

            segments = create_segments(chromosome_lengths[chrom], segment_length)

            segments['count'] = remixt.segalg.contained_counts(
                segments[['start', 'end']].values,
                reads[['start', 'end']].values,
            )

            segments['chromosome'] = chrom
            segments['num_obs'] = 1

            segments.to_csv(cna, sep='\t', index=False, header=False,
                columns=['chromosome', 'end', 'count', 'num_obs'])


def write_tumour_baf(baf_filename, normal_filename, tumour_filename):

    with open(baf_filename, 'w') as baf_file:

        chromosomes = remixt.seqdataio.read_chromosomes(normal_filename)

        for chrom in chromosomes:

            normal_allele_count = remixt.analysis.haplotype.read_snp_counts(normal_filename, chrom)

            remixt.analysis.haplotype.infer_snp_genotype(normal_allele_count)

            het_positions = normal_allele_count.loc[normal_allele_count['AB'] == 1, ['position']]

            tumour_allele_count = remixt.analysis.haplotype.read_snp_counts(tumour_filename, chrom)
            tumour_allele_count = tumour_allele_count.merge(het_positions)

            tumour_allele_count['ref_count'] = tumour_allele_count['ref_count'].astype(int)
            tumour_allele_count['alt_count'] = tumour_allele_count['alt_count'].astype(int)

            tumour_allele_count['minor_count'] = np.minimum(
                tumour_allele_count['ref_count'],
                tumour_allele_count['alt_count'],
            )

            tumour_allele_count['total_count'] = (
                tumour_allele_count['ref_count'] +
                tumour_allele_count['alt_count']
            )

            tumour_allele_count['chromosome'] = chrom

            tumour_allele_count.to_csv(baf_file, sep='\t', index=False, header=False,
                columns=['chromosome', 'position', 'minor_count', 'total_count'])


def _get_segments(filename):
    segments = set()
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split()
            segments.add((fields[0], fields[1]))
    return segments


def _filter_segments(filename, out_filename, segments):
    with open(filename) as f_in, open(out_filename, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
                continue
            fields = line.rstrip().split()
            if (fields[0], fields[1]) not in segments:
                continue
            f_out.write(line)


def _intersect_filtered(in_filenames, out_filenames):
    segments = None
    for in_filename in in_filenames:
        if segments is None:
            segments = _get_segments(in_filename)
        else:
            segments = segments.intersection(_get_segments(in_filename))

    for in_filename, out_filename in zip(in_filenames, out_filenames):
        _filter_segments(in_filename, out_filename, segments)


def prepare_data(
    normal_filename,
    tumour_filename,
    normal_cna_filename,
    tumour_cna_filename,
    tumour_baf_filename,
    config
):
    """ Initialize analysis
    """

    chromosome_lengths = read_chromosome_lengths(config['chrom_info_filename'])

    write_cna(normal_cna_filename, normal_filename, chromosome_lengths)
    write_cna(tumour_cna_filename, tumour_filename, chromosome_lengths)

    write_tumour_baf(tumour_baf_filename, normal_filename, tumour_filename)


def run_clonehd(
    normal_cna_filename,
    tumour_cna_filename,
    tumour_baf_filename,
    tumour_summary_filename,
    cna_subclone_filenames,
    baf_subclone_filenames,
    raw_data_dir,
):
    """ Run the analysis with specific initialization parameters

    """

    make_directory(raw_data_dir)

    num_clones = 2

    pypeliner.commandline.execute(
        'pre-filter',
        '--data', normal_cna_filename,
        '--pre', os.path.join(raw_data_dir, 'normal.cna'),
    )

    pypeliner.commandline.execute(
        'pre-filter',
        '--data', tumour_cna_filename,
        '--pre', os.path.join(raw_data_dir, 'tumour.cna'),
    )

    _intersect_filtered(
        (
            os.path.join(raw_data_dir, 'normal.cna.pref.txt'),
            os.path.join(raw_data_dir, 'tumour.cna.pref.txt'),
            tumour_baf_filename,
        ),
        (
            os.path.join(raw_data_dir, 'normal.cna.pref.intersect.txt'),
            os.path.join(raw_data_dir, 'tumour.cna.pref.intersect.txt'),
            os.path.join(raw_data_dir, 'tumour.baf.intersect.txt'),
        ),            
    )

    pypeliner.commandline.execute(
        'filterHD',
        '--data', os.path.join(raw_data_dir, 'normal.cna.pref.intersect.txt'),
        '--mode', '4',
        '--pre', os.path.join(raw_data_dir, 'normal.cna'),
        '--rnd', '0',
    )

    pypeliner.commandline.execute(
        'filterHD',
        '--data', os.path.join(raw_data_dir, 'tumour.cna.pref.intersect.txt'),
        '--mode', '4',
        '--pre', os.path.join(raw_data_dir, 'tumour.cna'),
        '--rnd', '0',
    )

    pypeliner.commandline.execute(
        'filterHD',
        '--data', os.path.join(raw_data_dir, 'tumour.cna.pref.intersect.txt'),
        '--mode', '4',
        '--pre', os.path.join(raw_data_dir, 'tumour.cna.bias'),
        '--bias', os.path.join(raw_data_dir, 'normal.cna.posterior-1.txt'),
        '--sigma', '0',
        '--jumps', '1',
        '--rnd', '0',
    )

    pypeliner.commandline.execute(
        'filterHD',
        '--data', os.path.join(raw_data_dir, 'tumour.baf.intersect.txt'),
        '--mode', '2',
        '--pre', os.path.join(raw_data_dir, 'tumour.baf'),
        '--sigma', '0',
        '--jumps', '1',
        '--reflect', '1',
        '--dist', '1',
        '--rnd', '0',
    )

    pypeliner.commandline.execute(
        'cloneHD',
        '--cna', os.path.join(raw_data_dir, 'tumour.cna.pref.intersect.txt'),
        '--baf', os.path.join(raw_data_dir, 'tumour.baf.intersect.txt'),
        '--pre', os.path.join(raw_data_dir, 'tumour'),
        '--bias', os.path.join(raw_data_dir, 'normal.cna.posterior-1.txt'),
        '--seed', '123',
        '--trials', '2',
        '--nmax', str(num_clones),
        '--force', str(num_clones),
        '--max-tcn', '6',
        '--cna-jumps', os.path.join(raw_data_dir, 'tumour.cna.bias.jumps.txt'),
        '--baf-jumps', os.path.join(raw_data_dir, 'tumour.baf.jumps.txt'),
        '--min-jump', '0.01',
        '--restarts', '10',
        '--mass-gauging', '1',
    )

    shutil.copyfile(os.path.join(raw_data_dir, 'tumour.summary.txt'), tumour_summary_filename)
    for clone_id in xrange(1, num_clones+1):
        shutil.copyfile(os.path.join(raw_data_dir, 'tumour.cna.subclone-{0}.txt'.format(clone_id)), cna_subclone_filenames[clone_id])
        shutil.copyfile(os.path.join(raw_data_dir, 'tumour.baf.subclone-{0}.txt'.format(clone_id)), baf_subclone_filenames[clone_id])


def report(
    tumour_summary_filename,
    cna_subclone_filenames,
    baf_subclone_filenames,
    results_filename,
    somatic_breakpoint_file=None,
):
    """ Report optimal copy number and mixture

    """

    segment_length = 1000

    with open(tumour_summary_filename, 'r') as summary_file:

        summary_info = dict()

        names = list()
        for line in summary_file:
            if line.startswith('#'):
                names = line[1:].split()
                if len(names) == 2 and names[1] == 'clones':
                    summary_info['num_clones'] = int(names[0])
                    names = ['mass'] + ['frac_'+str(i+1) for i in xrange(int(names[0]))]
            else:
                values = line.split()
                summary_info.update(dict(zip(names, values)))

    mix = [float(summary_info['frac_'+str(i+1)]) for i in xrange(summary_info['num_clones'])]
    mix = [1-sum(mix)] + mix

    cn_table = None

    for clone_id in xrange(1, summary_info['num_clones']+1):

        cna_filename = cna_subclone_filenames[clone_id]

        cna_data = pd.read_csv(cna_filename, delim_whitespace=True, converters={'#chr': str})
        cna_data.rename(columns={'#chr': 'chromosome', 'first-locus': 'start', 'last-locus': 'end'}, inplace=True)
        cna_data.drop(['nloci'], axis=1, inplace=True)
        # Unclear from the documentation, though techically the first datapoint is an endpoint
        # however, expanding regions results in inconsistencies
        cna_data['start'] -= segment_length
        cna_data.set_index(['chromosome', 'start', 'end'], inplace=True)
        cna_data = cna_data.idxmax(axis=1).astype(int)
        cna_data.name = 'total'
        cna_data = cna_data.reset_index()

        baf_filename = baf_subclone_filenames[clone_id]

        baf_data = pd.read_csv(baf_filename, delim_whitespace=True, converters={'#chr': str})
        baf_data.rename(columns={'#chr': 'chromosome', 'first-locus': 'start', 'last-locus': 'end'}, inplace=True)
        baf_data.drop(['nloci'], axis=1, inplace=True)
        baf_data.set_index(['chromosome', 'start', 'end'], inplace=True)
        baf_data = baf_data.fillna(0).idxmax(axis=1).astype(int)
        baf_data.name = 'allele'
        baf_data = baf_data.reset_index()

        data = remixt.segalg.reindex_segments(cna_data, baf_data)
        data = data.merge(cna_data[['total']], left_on='idx_1', right_index=True)
        data = data.merge(baf_data[['allele']], left_on='idx_2', right_index=True)

        data['major'] = np.maximum(data['allele'], data['total'] - data['allele'])
        data['minor'] = np.minimum(data['allele'], data['total'] - data['allele'])
        data.drop(['idx_1', 'idx_2', 'allele'], axis=1, inplace=True)

        # Having minor copies < 0 is common enough in the results
        # that we need to correct for it
        data['minor'] = np.maximum(data['minor'], 0)

        if not (data['minor'] >= 0).all():
            error = 'Negative minor copies\n'
            error += data[data['minor'] < 0].to_string()
            raise Exception(error)

        data.rename(inplace=True, columns={
            'total': 'total_{0}'.format(clone_id),
            'minor': 'minor_{0}'.format(clone_id),
            'major': 'major_{0}'.format(clone_id),
        })

        if cn_table is None:
            cn_table = data

        else:
            cn_table_prev = cn_table
            cn_table = remixt.segalg.reindex_segments(cn_table_prev, data)

            cn_table_prev.drop(['chromosome', 'start', 'end'], axis=1, inplace=True)
            data.drop(['chromosome', 'start', 'end'], axis=1, inplace=True)

            cn_table = cn_table.merge(cn_table_prev, left_on='idx_1', right_index=True)
            cn_table = cn_table.merge(data, left_on='idx_2', right_index=True)

            cn_table.drop(['idx_1', 'idx_2'], axis=1, inplace=True)

    # Post-hoc breakpoint copy number
    if somatic_breakpoint_file is not None:
        brk_cn = calculate_breakpoint_copy_number(somatic_breakpoint_file, cn_table)
    else:
        brk_cn = None

    with pd.HDFStore(results_filename, 'w') as store:
        store['mix'] = pd.Series(mix)
        store['cn'] = cn_table
        if brk_cn is not None:
            store['brk_cn'] = brk_cn



