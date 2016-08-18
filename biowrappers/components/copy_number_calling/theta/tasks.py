import os
import pandas as pd

import pypeliner.commandline
import remixt.seqdataio
import remixt.analysis.haplotype

# wget http://compbio.med.harvard.edu/BIC-seq/hg19.CRG.75bp.tar.gz
# wget ftp://ftp.ensembl.org/pub/release-70/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.70.dna.chromosome.20.fa.gz
# wget ftp://ftp.ensembl.org/pub/release-70/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.70.dna.chromosome.21.fa.gz

config = {}
config['chromosomes'] = ['20', '21']
config['chromosome_fasta'] = {
    '20': '/genesis/shahlab/amcpherson/theta_test/Homo_sapiens.GRCh37.70.dna.chromosome.20.fa',
    '21': '/genesis/shahlab/amcpherson/theta_test/Homo_sapiens.GRCh37.70.dna.chromosome.21.fa',
}
config['mappability_file'] = {
    '20': '/genesis/shahlab/amcpherson/theta_test/hg19.CRG.75bp/hg19.CRC.75mer.chr20.txt',
    '21': '/genesis/shahlab/amcpherson/theta_test/hg19.CRG.75bp/hg19.CRC.75mer.chr21.txt',
}

seg_output_filename = os.path.join(tmp_directory, 'seg.output')


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

    allele_count.to_csv(allele_filename, sep='\t', index=False, header=False)

## prefix for filenames TODO
def run_bicseq2_norm(prefix, seqdata_filename, config, tmp_directory, read_length, fragment_length):
    reads_filenames = {}
    norm_filenames = {}
    for chromosome in config['chromosomes']:
        reads_filenames[chromosome] = os.path.join(tmp_directory, 'chr{}.seq'.format(chromosome))
        norm_filenames[chromosome] = os.path.join(tmp_directory, 'chr{}.norm.bin'.format(chromosome))

    for chromosome in config['chromosomes']:
        chrom_reads = remixt.seqdataio.read_filtered_fragment_data(seqdata_filename, chromosome=chromosome)
        chrom_reads[['start']].to_csv(reads_template.format(chromosome), index=False, header=False)

    norm_config_filename = os.path.join(tmp_directory, 'norm.config')

    with open(norm_config_filename, 'w') as f:
        f.write('chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n')
        for chromosome in config['chromosomes']:
            f.write('{}\t'.format(chromosome))
            f.write('{}\t'.format(config['chromosome_fasta'][chromosome]))
            f.write('{}\t'.format(config['mappability_file'][chromosome]))
            f.write('{}\t'.format(reads_template.format(chromosome)))
            f.write('{}\n'.format(norm_filenames[chromosome]))

    pypeliner.commandline.execute(
        'bicseq2-norm',
        '-l', read_length,
        '-s', fragment_length,
        '-b', str(config['bin_size']),
        '--tmp', tmp_directory,
        norm_config_filename,
        norm_output_filename,
    )


def run_bicseq2_seg(seg_output_filename, normal_filename, tumour_filename, config, tmp_directory, read_length, fragment_length):
    normal_norm_filename = os.path.join(tmp_directory, 'normal.norm.output')
    tumour_norm_filename = os.path.join(tmp_directory, 'tumour.norm.output')

    normal_reads_template = os.path.join(tmp_directory, 'normal.chr{}.seq'.format(chromosome))
    normal_norm_template = os.path.join(tmp_directory, 'normal.chr{}.norm.bin'.format(chromosome))

    tumour_reads_template = os.path.join(tmp_directory, 'tumour.chr{}.seq'.format(chromosome))
    tumour_norm_template = os.path.join(tmp_directory, 'tumour.chr{}.norm.bin'.format(chromosome))

    run_bicseq2_norm(normal_norm_filename, normal_filename, config, tmp_directory, read_length, fragment_length)
    run_bicseq2_norm(tumour_norm_filename, tumour_filename, config, tmp_directory, read_length, fragment_length)

    seg_config_filename = os.path.join(tmp_directory, 'seg.config')

    with open(seg_config_filename, 'w') as f:
        f.write('chromName\binFileNorm.Case\tbinFileNorm.Control\n')
        for chromosome in config['chromosomes']:
            f.write('{}\t'.format(chromosome))
            f.write('{}\n'.format(tumour_norm_filenames[chromosome]))

    pypeliner.commandline.execute(
        'bicseq2-seg',
        '--control',
        '--tmp', tmp_directory,
        seg_config_filename,
        seg_output_filename,
    )


def run_theta(normal_filename, tumour_filename, config, tmp_directory, read_length, fragment_length):
    normal_allele_filename = os.path.join(tmp_directory, 'normal_alleles.tsv')
    tumour_allele_filename = os.path.join(tmp_directory, 'tumour_alleles.tsv')

    normal_allele_count = calculate_allele_counts(normal_filename, chromosomes=config['chromosomes'])
    tumour_allele_count = calculate_allele_counts(tumour_filename, chromosomes=config['chromosomes'])

    write_theta_format_alleles(normal_allele_filename, normal_allele_count)
    write_theta_format_alleles(tumour_allele_filename, tumour_allele_count)

    bicseq2_seg_filename = os.path.join(tmp_directory, 'bicseq2.seg')
    run_bicseq2_seg(bicseq2_seg_filename, normal_filename, tumour_filename, config, tmp_directory, read_length, fragment_length)

    theta_seg_filename = os.path.join(tmp_directory, 'theta.seg')
    pypeliner.commandline.execute(
        'BICSeqToTHetA',
        bicseq2_seg_filename,
        '-OUTPUT_PREFIX', theta_seg_filename,
    )

    theta_results_filename = os.path.join(tmp_directory, 'theta.txt')
    pypeliner.commandline.execute(
        'RunTHetA', '--FORCE',
        os.path.abspath(theta_seg_filename),
        '--TUMOR_FILE', os.path.abspath(normal_allele_filename),
        '--NORMAL_FILE', os.path.abspath(tumour_allele_filename),
        '--RESULTS', os.path.abspath(theta_results_filename),
    )


