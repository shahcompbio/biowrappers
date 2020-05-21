'''
Created on Nov 2, 2015

@author: Andrew Roth
'''
import biowrappers.components.variant_calling.utils as utils
import pandas as pd
import vcf
from bx.bbi.bigwig_file import BigWigFile


def get_mappability(
        mappability_file,
        vcf_file,
        out_file,
        region=None,
        append_chr=True):
    map_reader = BigWigFile(open(mappability_file, 'rb'))

    vcf_reader = vcf.Reader(filename=vcf_file)

    if region is not None:
        chrom, beg, end = utils.parse_region_for_vcf(region)
        try:
            vcf_reader = vcf_reader.fetch(chrom, start=beg, end=end)
        except ValueError:
            print("no data for region {} in vcf".format(region))
            vcf_reader = []

    data = []

    for record in vcf_reader:
        if append_chr:
            chrom = 'chr{0}'.format(record.CHROM)

        else:
            chrom = record.CHROM

        coord = record.POS

        beg = coord - 100

        beg = max(beg, 0)

        end = coord + 100

        result = map_reader.query(chrom, beg, end, 1)

        if result is None:
            mappability = 0

        else:
            mappability = result[0]['mean']

        data.append({'chrom': record.CHROM, 'coord': record.POS, 'mappability': mappability})

    data = pd.DataFrame(data)

    if out_file.endswith('.gz.tmp'):
        data.to_csv(out_file, index=False, compression='gzip')
    else:
        data.to_csv(out_file, index=False)
