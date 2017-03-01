'''
Created on Dec 13, 2015

@author: Andrew Roth
'''
import os
import pypeliner.app


def add_pypeliner_args(parser):
    pypeliner.app.add_arguments(parser)


def add_variant_calling_region_args(parser):
    from biowrappers.components.variant_calling.utils import default_chromosomes

    parser.add_argument('--chromosomes', nargs='+', default=default_chromosomes)

    parser.add_argument('--split_size', default=int(1e7), type=int)


def add_normal_tumour_bam_args(parser):
    parser.add_argument('-nb', '--normal_bam_file', required=True)

    parser.add_argument('-tb', '--tumour_bam_file', required=True)


def add_normal_multiple_tumour_bam_args(parser):
    parser.add_argument('-nb', '--normal_bam_file', required=True)

    parser.add_argument('-tb', '--tumour_bam_files', nargs='+', required=True)


def add_ref_genome_arg(parser):
    parser.add_argument('-rg', '--ref_genome_fasta_file', required=True)


def load_pypeliner_config(args):
    return vars(args)


def get_tumour_bam_file_dict(args):
    tumour_bam_files = {}
    for bam_file in args.tumour_bam_files:
        sample_id = os.path.basename(bam_file).rstrip('.bam')
        tumour_bam_files[sample_id] = bam_file

    return tumour_bam_files
