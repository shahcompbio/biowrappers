'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
import numpy as np
import pandas as pd
import pysam
import vcf

import biowrappers.components.variant_calling.utils as utils

nucleotides = ('A', 'C', 'G', 'T')

#=======================================================================================================================
# Allele counting
#=======================================================================================================================


def get_snv_allele_counts_for_vcf_targets(
        bam_file,
        vcf_file,
        out_file,
        table_name,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        region=None,
        vcf_to_bam_chrom_map=None):

    bam = pysam.AlignmentFile(bam_file, 'rb')

    vcf_reader = vcf.Reader(filename=vcf_file)

    if region is not None:
        chrom, beg, end = utils.parse_region_for_vcf(region)

        try:
            vcf_reader = vcf_reader.fetch(chrom, start=beg, end=end)

        except ValueError:
            vcf_reader = ()

    data = []

    for record in vcf_reader:
        if vcf_to_bam_chrom_map is not None:
            bam_chrom = vcf_to_bam_chrom_map[record.CHROM]

        else:
            bam_chrom = record.CHROM

        df = _get_counts_df(
            bam,
            bam_chrom,
            record.POS,
            record.POS + 1,
            count_duplicates=count_duplicates,
            min_bqual=min_bqual,
            min_mqual=min_mqual,
            strand='both'
        )

        counts = df.iloc[0]

        ref_base = record.REF

        # Skip record with reference base == N
        if ref_base not in nucleotides:
            continue

        for alt_base in record.ALT:
            alt_base = str(alt_base)

            if (len(ref_base) != 1) or (len(alt_base) != 1):
                continue

            # Skip record with alt base == N
            if alt_base not in nucleotides:
                continue

            # Format output record
            out_row = {
                'chrom': record.CHROM,
                'coord': record.POS,
                'ref': ref_base,
                'alt': alt_base,
                'ref_counts': counts[ref_base],
                'alt_counts': counts[alt_base]
            }

            data.append(out_row)

    data = pd.DataFrame(data, columns=['chrom', 'coord', 'ref', 'alt', 'ref_counts', 'alt_counts'])

    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    hdf_store[table_name] = data

    hdf_store.close()


def get_snv_allele_counts_for_region(
        bam_file,
        out_file,
        region,
        table_name,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        report_non_variant_positions=True,
        report_zero_count_positions=False):

    bam = pysam.AlignmentFile(bam_file, 'rb')

    chrom, beg, end = _parse_region(region)

    forward_df = _get_counts_df(
        bam,
        chrom,
        beg,
        end,
        count_duplicates=count_duplicates,
        min_bqual=min_bqual,
        min_mqual=min_mqual,
        strand='forward'
    )

    reverse_df = _get_counts_df(
        bam,
        chrom,
        beg,
        end,
        count_duplicates=count_duplicates,
        min_bqual=min_bqual,
        min_mqual=min_mqual,
        strand='reverse'
    )

    reverse_df = reverse_df.rename(columns=lambda x: x.lower())

    df = pd.concat([forward_df, reverse_df], axis=1)

    if not report_zero_count_positions:
        df = df[df.sum(axis=1) > 0]

    if (not report_non_variant_positions) and (df.shape[0] > 0):
        df = df[df.apply(_get_variant_positions_strand, axis=1)]

    df.reset_index(inplace=True)

    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    hdf_store[table_name] = df

    hdf_store.close()


def get_variant_position_counts(
        normal_bam_file,
        tumour_bam_files,
        out_file,
        region,
        table_group,
        count_duplicates=False,
        min_bqual=0,
        min_mqual=0,
        min_normal_depth=0,
        min_tumour_depth=0,
        min_variant_depth=0,
        report_strand_counts=True):
    """ Get counts for positions with at least two alleles in one or more tumour samples.

    This function filters for all positions which exceed the minimum depth in the normal sample and at least one tumour
    sample. It further filters for positions which have at least two alleles present in one or more tumour samples.

    """

    chrom, beg, end = _parse_region(region)

    tumour_samples = tumour_bam_files.keys()

    bams = {'normal': pysam.AlignmentFile(normal_bam_file, 'rb')}

    for sample in tumour_samples:
        bams[sample] = pysam.AlignmentFile(tumour_bam_files[sample], 'rb')

    counts = {}

    for sample in bams:
        if report_strand_counts:
            forward_df = _get_counts_df(
                bams[sample],
                chrom,
                beg,
                end,
                count_duplicates=count_duplicates,
                min_bqual=min_bqual,
                min_mqual=min_mqual,
                strand='forward'
            )

            reverse_df = _get_counts_df(
                bams[sample],
                chrom,
                beg,
                end,
                count_duplicates=count_duplicates,
                min_bqual=min_bqual,
                min_mqual=min_mqual,
                strand='reverse'
            )

            reverse_df.rename(columns=lambda x: x.lower(), inplace=True)

            counts[sample] = pd.concat([forward_df, reverse_df], axis=1)

        else:
            counts[sample] = _get_counts_df(
                bams[sample],
                chrom,
                beg,
                end,
                count_duplicates=count_duplicates,
                min_bqual=min_bqual,
                min_mqual=min_mqual,
                strand='both'
            )

    # Depth filtering
    valid_positions = np.zeros(counts['normal'].shape[0], dtype=bool)

    for sample in tumour_samples:
        valid_positions = np.logical_or(valid_positions, counts[sample].sum(axis=1) >= min_tumour_depth)

    valid_positions = np.logical_and(valid_positions, counts['normal'].sum(axis=1) >= min_normal_depth)

    for sample in counts:
        counts[sample] = counts[sample][valid_positions]

    # Skip filtering if there are no valid positions
    if np.sum(valid_positions) > 0:
        # Variant position filtering
        variant_positions = np.zeros(counts['normal'].shape[0], dtype=bool)

        for sample in tumour_samples:
            if report_strand_counts:
                variant_positions = np.logical_or(
                    variant_positions, counts[sample].apply(_get_variant_positions_strand, axis=1))

            else:
                variant_positions = np.logical_or(
                    variant_positions, counts[sample].apply(lambda x: _get_variant_positions(x, min_variant_depth), axis=1))

        for sample in counts:
            counts[sample] = counts[sample][variant_positions]

    # Write output
    hdf_store = pd.HDFStore(out_file, 'w', complevel=9, complib='blosc')

    for sample in counts:
        df = counts[sample]

        df.reset_index(inplace=True)

        table_name = '/'.join((table_group, sample))

        hdf_store.append(table_name, df)

    hdf_store.close()


def _get_counts_df(bam_file,
                   chrom,
                   beg,
                   end,
                   count_duplicates=False,
                   min_bqual=30,
                   min_mqual=30,
                   strand='both'):
    '''
    Get counts 1 based indexing.
    '''

    x = bam_file.count_coverage(
        chrom,
        beg - 1,
        end - 1,
        quality_threshold=min_bqual,
        read_callback=lambda x: _check_read(
            x,
            count_duplicates=count_duplicates,
            min_mqual=min_mqual,
            strand=strand)
    )

    df = pd.DataFrame(pd.np.array(x).T, columns=['A', 'C', 'G', 'T'])

    df.insert(0, 'chrom', chrom)

    df.insert(1, 'coord', np.arange(beg, end))

    df = df.set_index(['chrom', 'coord'])

    return df


def _check_read(read, count_duplicates=False, min_mqual=30, strand='both'):
    valid = True

    if read.mapping_quality < min_mqual:
        valid = False

    elif read.is_duplicate and (not count_duplicates):
        valid = False

    elif read.is_unmapped:
        valid = False

    elif read.is_qcfail:
        valid = False

    elif read.is_secondary:
        valid = False

    elif (strand == 'reverse') and (not read.is_reverse):
        valid = False

    elif (strand == 'forward') and (read.is_reverse):
        valid = False

    return valid


def _get_variant_positions(row, min_variant_depth):
    counts = sorted(row, reverse=True)

    if counts[1] > min_variant_depth:
        return True

    else:
        return False


def _get_variant_positions_strand(row):
    counts = [row[x.lower()] + row[x.upper()] for x in nucleotides]

    counts = sorted(counts, reverse=True)

    if counts[1] > 0:
        return True

    else:
        return False


def _parse_region(region):
    chrom, coords = region.split(':')

    beg, end = coords.split('-')

    beg = int(beg)

    end = int(end)

    return chrom, beg, end
