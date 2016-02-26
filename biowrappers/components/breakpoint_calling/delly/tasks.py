import vcf
import pandas as pd


def convert_vcf(vcf_filename, store_filename):
    vcf_reader = vcf.Reader(filename=vcf_filename)

    breakpoint_table = list()
    breakpoint_library_table = list()

    for row in vcf_reader:
        prediction_id = row.ID

        chrom_1 = row.CHROM
        chrom_2 = row.INFO['CHR2']

        strand_1, strand_2 = [('-', '+')[a == '3'] for a in row.INFO['CT'].split('to')]

        coord_1 = row.POS
        coord_2 = row.sv_end

        if 'LowQual' in row.FILTER:
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

