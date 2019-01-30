import pandas as pd
import remixt.seqdataio
import remixt.analysis.haplotype


def calculate_allele_counts(seqdata_filename):
    """ Calculate allele counts from seqdata.
    """
    allele_counts = [pd.DataFrame(columns = ["position", "ref_count", "alt_count", "chromosome"])]

    chromosomes = remixt.seqdataio.read_chromosomes(seqdata_filename)

    for chrom in chromosomes:
        chrom_allele_counts = remixt.analysis.haplotype.read_snp_counts(seqdata_filename, chrom)
        chrom_allele_counts['chromosome'] = chrom
        allele_counts.append(chrom_allele_counts)

    if len(allele_counts) > 0:
        allele_counts = pd.concat(allele_counts, ignore_index=True)
    else:
        allele_counts = pd.DataFrame()

    allele_counts['ref_count'] = allele_counts['ref_count'].astype(int)
    allele_counts['alt_count'] = allele_counts['alt_count'].astype(int)
    allele_counts['position'] = allele_counts['position'].astype(int)

    return allele_counts


