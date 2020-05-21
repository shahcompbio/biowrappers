'''
Created on Nov 2, 2015

@author: Andrew Roth
'''
import pandas as pd
import vcf


def annotate_db_status(db_vcf_file, target_vcf_file, out_file):
    db_reader = vcf.Reader(filename=db_vcf_file)

    reader = vcf.Reader(filename=target_vcf_file)

    data = []

    for record in reader:
        chrom = record.CHROM

        coord = record.POS

        try:
            db_position_records = [x for x in db_reader.fetch(chrom, coord - 1, coord)]

        except ValueError:
            db_position_records = []

        for db_record in db_position_records:

            if (db_record.CHROM != chrom) or (db_record.POS != coord):
                continue

            if db_record.is_indel:
                indel = 1

            else:
                indel = 0

            for alt in record.ALT:

                if (record.REF == db_record.REF) and (alt in db_record.ALT):
                    exact_match = 1

                else:
                    exact_match = 0

                out_row = {
                    'chrom': chrom,
                    'coord': coord,
                    'ref': record.REF,
                    'alt': str(alt),
                    'db_id': db_record.ID,
                    'exact_match': exact_match,
                    'indel': indel
                }

                data.append(out_row)

    data = pd.DataFrame(data)

    if out_file.endswith('.gz.tmp'):
        data.to_csv(out_file, index=False, compression='gzip')
    else:
        data.to_csv(out_file, index=False)
