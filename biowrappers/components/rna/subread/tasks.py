import pypeliner


def run_feature_counts(in_file, gene_gtf_file,  out_file):
    cmd = [
        'featureCounts',
        '-a', gene_gtf_file,
        '-o', out_file,
        in_file
    ]
    pypeliner.commandline.execute(*cmd)
