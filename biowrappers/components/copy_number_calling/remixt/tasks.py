import pandas as pd


def select_solution(selected_file, results_file):
    common_tables = [
        '/breakpoint_adjacencies',
        '/breakpoints',
        '/minor_modes',
        '/read_depth',
        '/reference_adjacencies',
    ]

    solution_tables = [
        '/betabin_M',
        '/brk_cn',
        '/cn',
        '/h',
        '/h_init',
        '/negbin_r',
    ]

    with pd.HDFStore(selected_file, 'w') as selected_store, pd.HDFStore(results_file, 'r') as results_store:
        for table_name in common_tables:
            selected_store[table_name] = results_store[table_name]

        stats = results_store['stats']
        solution_stats = stats.loc[stats['bic_optimal']].iloc[0]
        solution_idx = solution_stats['idx']

        for table_name in solution_tables:
            selected_store[table_name] = results_store['/solutions/solution_{}/{}'.format(solution_idx, table_name)]

