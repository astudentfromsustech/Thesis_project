import pathlib as p
import pandas as pd
import numpy as np

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    nums = ['WT', 'HS']
    for num in nums:
        print(num)
        src_dir = root.parent / 'results' / 'PNAS_data_V2_results' / f'{num}_combined_10000_targetChroms_normalize0_1_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'

        TAD_raw = pd.read_csv(src_dir / f'{num}_domains.bed', header=None, sep='\t')
        print(TAD_raw.shape)
        print(TAD_raw.head())
        df1 = TAD_raw.iloc[:, 0:3]
        TAD = pd.concat([df1, df1], axis=1)
        print(TAD.shape)
        print(TAD.head())
        TAD.to_csv(src_dir /  f'{num}_domains.bedpe', header=None, index=False, sep='\t')