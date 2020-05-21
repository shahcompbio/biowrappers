import pandas as pd


def concatenate_csv(in_filenames, out_filename):
    if isinstance(in_filenames, dict):
        in_filenames = in_filenames.values()

    data = []

    for in_filename in in_filenames:
        try:
            df = pd.read_csv(in_filename)
            data.append(df)
        except pd.errors.EmptyDataError:
            pass

    data = pd.concat(data)

    if out_filename.endswith('.gz.tmp'):
        data.to_csv(out_filename, index=False, compression='gzip')
    else:
        data.to_csv(out_filename, index=False)

