import pathlib as p
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / '1_pairs_from_dm6_to_dm3_usingChunks_combineChunks'
    anno_dir = root.parent / 'annotations'
    dest_dir = root.parent / 'data' / '2_pairs_downsample'
    dest_dir.mkdir(exist_ok=True, parents=True)

    chrom_order = ['2L', '2R', '3L', '3R', '4', 'X']
    pair_nums = [10000000, 20000000, 40000000, 60000000, 80000000, 100000000]
    conditions = ['WT', 'HS']
    # conditions = ['HS']

    for condition in conditions:
        print(condition)
        data = pd.read_csv(src_dir / f'{condition}_combined_PNAS_chunksCombined.pairs.Pre', header=None, sep=' ')
        data.columns = ['strand1', 'chrom1_dm3', 'pos1_dm3', 'frag1', 'strand2', 'chrom2_dm3', 'pos2_dm3', 'frag2']
        data['chrom1_dm3'] = pd.Categorical(data['chrom1_dm3'], categories=chrom_order, ordered=True)
        data['chrom2_dm3'] = pd.Categorical(data['chrom2_dm3'], categories=chrom_order, ordered=True)
        data['pos1_dm3'] = data['pos1_dm3'].astype(int)
        data['pos2_dm3'] = data['pos2_dm3'].astype(int)
        print(data.shape)
        for pair_num in pair_nums:
            print(pair_num)
            data_dn = data.sample(n=pair_num, replace=False)
            data_dn_sorted = data_dn.sort_values(by=['chrom1_dm3', 'chrom2_dm3'])
            print(data_dn_sorted.shape)
            print(data_dn_sorted.head())
            data_dn_sorted.to_csv(dest_dir / f'{condition}_combined_PNAS.pairs.downsample{pair_num}.Pre', index=None, header=None, sep=' ')




