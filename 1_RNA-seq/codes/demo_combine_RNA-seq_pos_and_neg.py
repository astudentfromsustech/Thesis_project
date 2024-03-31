import pathlib as p
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)



if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data'
    dest_dir = root.parent / 'results'
    dest_dir.mkdir(parents=True, exist_ok=True)

    num = 'RDS029'
    neg = pd.read_csv(src_dir / f'{num}NEW.negative.bedGraph', sep='\t', header=None)
    neg.columns = ['chromosome', 'start', 'end', 'count']
    neg['count'] = neg['count'] * (-1)
    print(neg.shape)
    print(neg.head())

    
    pos = pd.read_csv(src_dir / f'{num}NEW.positive.bedGraph', sep='\t', header=None)
    pos.columns = ['chromosome', 'start', 'end', 'count']
    print(pos.shape)
    print(pos.head())

    stacked_df = pd.concat([pos, neg]).sort_values(by=['chromosome', 'start'])
    print(stacked_df.shape)
    print(stacked_df.head())
    stacked_df.to_csv(src_dir / f'{num}NEW.pos.neg.bedGraph', header=None, index=None, sep='\t')







