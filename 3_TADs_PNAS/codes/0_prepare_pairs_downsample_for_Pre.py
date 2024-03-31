import pathlib as p
import subprocess
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'pairs_downsample'
    dest_dir = root.parent / 'data' / '0_prepare_pairs_downsample_for_Pre'
    dest_dir.mkdir(exist_ok=True, parents=True)

    conditions = ['WT', 'HS']
    # conditions = ['WT']
    # pair_nums= [10000000, 20000000, 40000000, 60000000, 80000000, 100000000]
    pair_nums = [60000000, 80000000, 100000000]

    for condition in conditions:
        print(condition)
        for pair_num in pair_nums:
            print(pair_num)

            raw_data = pd.read_csv(src_dir / f'{condition}_combined_PNAS.pairs_dn{pair_num}', header=None, sep='\t')
            raw_data.columns = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type', 'mapq1', 'mapq2']
            raw_data['strand1'] = raw_data['strand1'].replace({'-': '0', '+': '1'})
            raw_data['strand2'] = raw_data['strand2'].replace({'-': '0', '+': '1'})
            # raw_data['chrom1'] = raw_data['chrom1'].str.replace('chr', '')
            # raw_data['chrom2'] = raw_data['chrom2'].str.replace('chr', '')
            print(raw_data.shape)
            print(raw_data.head())
            data = raw_data[['strand1', 'chrom1', 'pos1', 'strand2', 'chrom2', 'pos2', ]]
            data.insert(3, 'frag1', 0)
            data.insert(7, 'frag2', 1)
            print(data.shape)
            print(data.head())
            data_sorted = data.sort_values(by=['chrom1', 'chrom2', 'pos1', 'pos2'])
            print(data_sorted.shape)
            print(data_sorted.head())
            data_sorted.to_csv(dest_dir / f'{condition}_combined_PNAS.pairs_dn{pair_num}.Pre', index=None, header=None, sep=' ')

