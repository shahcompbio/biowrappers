import os
import pysam
import pandas as pd

import pypeliner.commandline


def write_samples_table(sample_types, sample_table_filename):
    with open(sample_table_filename, 'w') as f:
        for sample_id, sample_type in sample_types.iteritems():
            f.write('{}\t{}\n'.format(sample_id, sample_type))


def _rename_index(out_file):
    if out_file.endswith('.tmp'):
        out_index = out_file + '.csi'
        renamed_out_index = out_file[:-4] + '.csi'
        try:
            os.remove(renamed_out_index)
        except OSError:
            pass
        os.rename(out_index, renamed_out_index)


def run_delly_call(sv_type, delly_excl_chrom, ref_genome_fasta_file, bam_files, out_file):
    delly_args = [
        'delly', 'call',
        '-t', sv_type,
        '-x', delly_excl_chrom,
        '-g', ref_genome_fasta_file,
        '-o', out_file,
    ]

    delly_args += bam_files
    
    pypeliner.commandline.execute(*delly_args)
    _rename_index(out_file)


def run_delly_filter(sv_type, sample_file, ref_genome_fasta_file, in_file, out_file):
    delly_args = [
        'delly', 'filter',
        '-t', sv_type,
        '-f', 'somatic',
        '-o', out_file,
        '-s', sample_file,
        '-g', ref_genome_fasta_file,
        in_file,
    ]

    pypeliner.commandline.execute(*delly_args)
    _rename_index(out_file)


def convert_vcf(bcf_filename, store_filename):
    bcf_reader = pysam.VariantFile(bcf_filename, 'rb')

    breakpoint_table = list()
    breakpoint_library_table = list()

    for row in bcf_reader:
        prediction_id = row.id

        chrom_1 = row.chrom
        chrom_2 = row.info['CHR2']

        strand_1, strand_2 = [('-', '+')[a == '3'] for a in row.info['CT'].split('to')]

        coord_1 = row.pos
        coord_2 = row.info['END']

        if 'LowQual' in row.filter:
            qual = 0
        else:
            qual = 1

        breakpoint_table.append((prediction_id, chrom_1, chrom_2, strand_1, strand_2, coord_1, coord_2, qual))

        for call in row.samples:
            library = call.sample

            num_spanning = call.data.DV
            num_split = call.data.RV

            breakpoint_library_table.append((prediction_id, library, num_spanning, num_split))

    breakpoint_table = pd.DataFrame(
        breakpoint_table,
        columns=[
            'prediction_id',
            'chromosome_1', 'chromosome_2',
            'strand_1', 'strand_2',
            'position_1', 'position_2',
            'qual',
        ]
    )

    breakpoint_library_table = pd.DataFrame(
        breakpoint_library_table,
        columns=[
            'prediction_id',
            'library',
            'num_spanning',
            'num_split'
        ]
    )

    with pd.HDFStore(store_filename, 'w') as store:
        store['/breakpoint'] = breakpoint_table
        store['/breakpoint_library'] = breakpoint_library_table

