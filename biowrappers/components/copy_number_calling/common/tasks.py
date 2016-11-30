import itertools
import numpy as np
import pandas as pd

import remixt.analysis.experiment
import remixt.analysis.pipeline
import remixt.cn_model


def calculate_breakpoint_copy_number(breakpoints_filename, cn_table, max_brk_dist=2000, max_seg_gap=int(3e6)):
    """ Post-hoc analysis of copy number transitions coincident with breakpoints.

    Args:
        breakpoints_filename (str): filename table of breakpoints in tsv format
        cn_table (pandas.DataFrame): table of segment copy number

    KwArgs:
        max_brk_dist (int): max distance from transition to breakend
        max_seg_gap (int): max gap between segments for wild type adjacencies

    """

    N = len(cn_table.index)

    if 'major_1' in cn_table:
        cn = [np.ones((N, 2), dtype=int)]
        for m in itertools.count(1):
            try:
                cn.append(cn_table[['major_{}'.format(m), 'minor_{}'.format(m)]].values)
            except KeyError:
                break
        cn = np.array(cn).swapaxes(0, 1)

    else:
        cn = [np.ones((N, 1), dtype=int) * 2]
        for m in itertools.count(1):
            try:
                cn.append(cn_table[['total_{}'.format(m)]].values)
            except KeyError:
                break
        cn = np.array(cn).swapaxes(0, 1)

    adjacencies = remixt.analysis.experiment.get_wild_type_adjacencies(cn_table, max_seg_gap)
    breakpoint_data = pd.read_csv(breakpoints_filename, sep='\t', converters={'chromosome_1': str, 'chromosome_2': str})
    breakpoint_segment_data = remixt.analysis.experiment.create_breakpoint_segment_table(cn_table, breakpoint_data, adjacencies, max_brk_dist=max_brk_dist)

    # Convert to a unique set of breakpoints
    breakpoints = set()
    for n_1, side_1, n_2, side_2 in breakpoint_segment_data[['n_1', 'side_1', 'n_2', 'side_2']].values:
        breakpoints.add(frozenset([(n_1, side_1), (n_2, side_2)]))

    breakpoints = list(breakpoints)

    brk_cn = remixt.cn_model.decode_breakpoints_naive(cn, adjacencies, breakpoints)
    brk_cn_table = remixt.analysis.experiment.create_brk_cn_table(brk_cn, breakpoint_segment_data)

    return brk_cn_table


