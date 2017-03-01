import pypeliner

from biowrappers.components.variant_calling.utils import default_chromosomes

import biowrappers.components.variant_calling.samtools as samtools


def main(args):
    native_spec = '-V -q all.q -l mem_token={mem}G,mem_free={mem}G,h_vmem={mem}G'

    config = {
        'tmpdir': args.log_dir,
        'pretend': False,
        'submit': 'asyncqsub',
        'nativespec': native_spec,
        'maxjobs': 100,
        'nocleanup': False
    }

    pyp = pypeliner.app.Pypeline([], config)

    scheduler = pyp.sch

    indel_vcf_file = args.out_prefix + '.indel.vcf.gz'

    snv_vcf_file = args.out_prefix + '.snv.vcf.gz'

    samtools.create_samtools_variant_calling_workflow(
        scheduler,
        args.bam_file,
        args.ref_genome_fasta_file,
        indel_vcf_file,
        snv_vcf_file,
        chromosomes=args.chromosomes,
        split_size=args.split_size
    )

    pyp.run()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--bam_file', required=True)

    parser.add_argument('--ref_genome_fasta_file', required=True)

    parser.add_argument('--out_prefix', required=True)

    parser.add_argument('--chromosomes', nargs='+', default=default_chromosomes)

    parser.add_argument('--log_dir', default='./')

    parser.add_argument('--split_size', default=int(1e7), type=int)

    args = parser.parse_args()

    main(args)
