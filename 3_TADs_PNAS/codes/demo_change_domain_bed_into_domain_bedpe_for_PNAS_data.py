import pathlib as p
import pandas as pd
import numpy as np

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    nums = ['WT', 'HS']
    res = 10000
    for num in nums:
        print(num)
        src_dir = root.parent / 'results' / f'{num}_combined_{res}_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.1_delta0.01_dfr'

        TAD_raw = pd.read_csv(src_dir / f'{num}_domains.bed', header=None, sep='\t')
        print(TAD_raw.shape)
        print(TAD_raw.head())
        df1 = TAD_raw.iloc[:, 0:3]
        TAD = pd.concat([df1, df1], axis=1)
        print(TAD.shape)
        print(TAD.head())
        TAD.to_csv(src_dir /  f'{num}_res{res}_domains.bedpe', header=None, index=False, sep='\t')