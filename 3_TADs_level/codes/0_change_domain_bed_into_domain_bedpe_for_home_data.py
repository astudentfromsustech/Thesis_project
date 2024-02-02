import pathlib as p
import pandas as pd
import numpy as np

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    nums = ['01', '02', '04', '05', '07', '08']
    thresh=0.1
    delta=0.01
    for num in nums:
        print(num)
        src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh{thresh}_delta{delta}_dfr'

        TAD_raw = pd.read_csv(src_dir / f'CDS0{num}D_domains.bed', header=None, sep='\t')
        print(TAD_raw.shape)
        print(TAD_raw.head())
        df1 = TAD_raw.iloc[:, 0:3]
        TAD = pd.concat([df1, df1], axis=1)
        print(TAD.shape)
        print(TAD.head())
        TAD.to_csv(src_dir /  f'CDS0{num}D_thresh{thresh}_delta{delta}_domains.bedpe', header=None, index=False, sep='\t')