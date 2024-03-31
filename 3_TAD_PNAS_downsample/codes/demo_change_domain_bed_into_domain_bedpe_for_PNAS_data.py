import pathlib as p
import pandas as pd
import numpy as np

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    conditions = ['WT', 'HS']
    pair_nums = [10000000, 20000000, 40000000, 60000000, 80000000, 100000000]
    res = 10000
    for condition in conditions:
        print(condition)
        for pair_num in pair_nums:
            print(pair_num)
            src_dir = root.parent / 'data' /f'9_TADs_resolution_{res}' / f'{condition}_combined_PNAS.pairs_dn{pair_num}_min30k_max100k_step10k_thresh0.1_delta0.01_dfr'

            TAD_raw = pd.read_csv(src_dir / f'{condition}_combined_PNAS.pairs_dn{pair_num}_domains.bed', header=None, sep='\t')
            print(TAD_raw.shape)
            print(TAD_raw.head())
            df1 = TAD_raw.iloc[:, 0:3]
            TAD = pd.concat([df1, df1], axis=1)
            print(TAD.shape)
            print(TAD.head())
            TAD.to_csv(src_dir /  f'{condition}_combined_PNAS.pairs_dn{pair_num}_domains.bedpe', header=None, index=False, sep='\t')