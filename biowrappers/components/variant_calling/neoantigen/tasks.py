import os
import shutil
import errno
import pypeliner.commandline


def remove_temp(temp_space):
    try:
        shutil.rmtree(temp_space)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

def run_optitype(reads_1_fastq, reads_2_fastq, hla_type_file, temp_space):
    remove_temp(temp_space)

    pypeliner.commandline.execute(
        'OptiTypePipeline.py', '-i', reads_1_fastq, reads_2_fastq,
        '--dna', '-v', '-o', temp_space)

    results_subdir = os.listdir(temp_space)

    assert len(results_subdir) == 1

    date_time_subdir = results_subdir[0]

    results_file = os.path.join(temp_space, date_time_subdir, date_time_subdir + '_result.tsv')
    os.rename(results_file, hla_type_file)


def run_pvacseq(vcf_file, hla_type_file, results_file, temp_space, config):
    remove_temp(temp_space)

    iedb_install_dir = config['iedb_install_dir']
    epitope_length = config['epitope_length']
    algorithms = config['binding_algorithms']

    hla_type = pd.read_csv(hla_type_file, sep='\t')
    hla_alleles = hla_type[['A1', 'A2', 'B1', 'B2', 'C1', 'C2']].iloc[0].values
    hla_alleles_arg = ','.join(['HLA-' + a for a in hla_alleles])

    results = {}

    for algorithm in algorithms:
        remove_temp(temp_space)

        pypeliner.commandline.execute(
            'pvacseq', 'run', vcf_file, 'tumour_sample', hla_alleles_arg, algorithm, 
            temp_space, '--iedb-install-directory', iedb_install_dir, '-e', epitope_length)

        results_filename = os.path.join(temp_space, 'MHC_Class_I', 'tumour_sample.final.tsv')
        results[algorithm] = pd.read_csv(results_filename, sep='\t')

    with pd.HDFStore(results_file, 'w') as results_store:
        for algorithm in algorithms:
            results_store[algorithm] = results[algorithm]

