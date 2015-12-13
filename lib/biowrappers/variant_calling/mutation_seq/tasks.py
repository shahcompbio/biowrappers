'''
Created on Oct 31, 2015

@author: Andrew Roth
'''
import csv
import gzip
import os
import pypeliner
import vcf

def run_mutation_seq(
    install_dir,
    normal_bam_file,
    tumour_bam_file,
    ref_genome_fasta_file,
    out_file,
    max_normal_variant_proportion=100,
    min_depth=0,
    min_tumour_variant_count=0,
    region=None,
    tumour_content=None,
    use_heuristic_filters=False,
    variant_probability_threshold=0):
    
    classify_script = os.path.join(install_dir, 'classify.py')
    
    model_file = os.path.join(install_dir, 'model_v4.1.2.npz')
    
    vcf_metadata_file = os.path.join(install_dir, 'metadata.config')

    cmd = [
        'python',
        classify_script,
        'normal:{0}'.format(normal_bam_file),
        'tumour:{0}'.format(tumour_bam_file),
        'reference:{0}'.format(ref_genome_fasta_file),
        'model:{0}'.format(model_file),
        '--config', vcf_metadata_file,
        '--coverage', min_depth,
        '--normal_variant', max_normal_variant_proportion,
        '--out', out_file,
        '--threshold', variant_probability_threshold,
        '--tumour_variant', min_tumour_variant_count
    ]
    
    if region is not None:
        cmd.extend(['--interval', region])
    
    if not use_heuristic_filters:
        cmd.append('--no_filter')
    
    if tumour_content is not None:
        cmd.extend(['--purity ', tumour_content])
    
    pypeliner.commandline.execute(*cmd)

def filter_snvs_by_threshold(in_file, out_file, threshold=0.85):
    
    reader = vcf.Reader(filename=in_file)
    
    with open(out_file, 'w') as out_fh:
        writer = vcf.Writer(out_fh, reader)
        
        for record in reader:
            if record.INFO['PR'] < threshold:
                continue
            
            if len(record.FILTER) > 0:
                continue
            
            writer.write_record(record)
            
def build_results_table(in_file, out_file):
    
    cols = (
        'chrom',
        'coord',
        'ref',
        'alt',
        'probability_somatic',
        'ref_counts_normal',
        'alt_counts_normal',
        'ref_counts_tumour',
        'alt_counts_tumour',
        'deletion_counts',
        'insertion_counts',
        'tri_nucleotide_context'
    )
    
    reader = vcf.Reader(filename=in_file)
    
    with gzip.GzipFile(out_file, 'w') as out_fh:
        writer = csv.DictWriter(out_fh, cols, delimiter='\t')
        
        writer.writeheader()
        
        for record in reader:
            for alt in record.ALT:
                row = {
                   'chrom' : record.CHROM,
                   'coord' : record.POS,
                   'ref' : record.REF,
                   'alt' : str(alt),
                   'probability_somatic' : record.INFO['PR'],
                   'ref_counts_normal' : record.INFO['NR'],
                   'alt_counts_normal' : record.INFO['NA'],
                   'ref_counts_tumour' : record.INFO['TR'],
                   'alt_counts_tumour' : record.INFO['TA'],
                   'deletion_counts' : record.INFO['ND'],
                   'insertion_counts' : record.INFO['NI'],
                   'tri_nucleotide_context' : record.INFO['TC']
                }
            
                writer.writerow(row)
