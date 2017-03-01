import pandas as pd


def extract_somatic_breakpoint(breakpoint_results, somatic_breakpoints_file, config):
    store = pd.HDFStore(breakpoint_results, 'r')

    breakpoints = store['/breakpoints/destruct/breakpoint']

    breakpoints = breakpoints[breakpoints['num_reads'] >= config['breakpoint_filtering']['minimum_num_reads']]

    breakpoints = breakpoints[[
        'prediction_id',
        'chromosome_1',
        'strand_1',
        'position_1',
        'chromosome_2',
        'strand_2',
        'position_2',
    ]]

    breakpoints.to_csv(somatic_breakpoints_file, sep='\t', index=False)
