import gzip
import pandas as pd

def main(args):
    cells = []
    
    for file_name in args.in_files:
        df = pd.read_csv(file_name, compression='gzip', sep='\t')
        
        cells.append(set(df['cell_id']))
    
    cells = list(set.intersection(*cells))

    data = []
    
    for cell_id in cells:
        data.append({'cell_id' : cell_id, 'sample' : cell_id.split(':')[0]})
    
    data = pd.DataFrame(data)
    
    with gzip.GzipFile(args.out_file, 'w') as fh:
        data.to_csv(fh, index=False, sep='\t')
    
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--in_files', nargs='+', required=True)
  
    parser.add_argument('--out_file', required=True)
    
    args = parser.parse_args()
    
    main(args)