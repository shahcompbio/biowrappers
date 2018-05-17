import pandas as pd
import remixt.seqdataio
import remixt.analysis.haplotype


def calculate_allele_counts(seqdata_filename):
    """ Calculate allele counts from seqdata.
    """
    allele_counts = list()

    chromosomes = remixt.seqdataio.read_chromosomes(seqdata_filename)

    for chrom in chromosomes:
        chrom_allele_counts = remixt.analysis.haplotype.read_snp_counts(seqdata_filename, chrom)
        chrom_allele_counts['chromosome'] = chrom
        allele_counts.append(chrom_allele_counts)

    if len(allele_counts) > 0:
        allele_counts = pd.concat(allele_counts, ignore_index=True)
    else:
        allele_counts = pd.DataFrame()

    return allele_counts


