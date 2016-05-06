import pandas as pd


def select_solution(selected_file, results_file, config):
    with pd.HDFStore(selected_file, 'w') as selected_store, pd.HDFStore(results_file, 'r') as results_store:
        stats = results_store['stats']
        stats = stats[stats['proportion_divergent'] < config['max_prop_diverge']].copy()
        stats.sort_values('log_likelihood', ascending=False, inplace=True)
        solution_idx = stats.loc[stats.index[0], 'init_id']

        for table_name in results_store.keys():
            if table_name.startswith('/solutions/solution_'):
                solution_table_prefix = '/solutions/solution_{}'.format(solution_idx)

                if table_name.startswith(solution_table_prefix):
                    sub_table_name = table_name[len(solution_table_prefix):]

                    selected_store[sub_table_name] = results_store[table_name]

            else:
                selected_store[table_name] = results_store[table_name]


