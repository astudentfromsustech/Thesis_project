import pathlib as p
import pandas as pd
import numpy as np

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'ac_me3_downsample_combined'
    dest_dir = root.parent / 'data' / 'ac_me3_downsample_combined_hic'
    dest_dir.mkdir(parents=True, exist_ok=True)

    WT_me3 = pd.read_csv(src_dir / f'CDS001D_combined.e500.juice', header=None, sep='\s+')
    print(WT_me3.shape)
    print(WT_me3.head())


    WT_ac = pd.read_csv(src_dir / f'CDS004D10DC.e500.juice', header=None, sep='\s+')
    print(WT_ac.shape)
    print(WT_ac.head())

    WT_me3_ac_combine = pd.concat([WT_me3, WT_ac]).sort_values(by=[1, 5])
    print(WT_me3_ac_combine.shape)
    print(WT_me3_ac_combine.head())

    # WT_me3_ac_combine.to_csv(src_dir / f'WT_me3_ac_combined.e500.juice', header=None, index=None, sep=' ')
    # WT_me3_ac_combine = pd.read_csv(src_dir / f'WT_me3_ac_combined.e500.juice', header=None, sep='\s+',
    #                                 dtype={1: str, 5: str})
    #
    # print(WT_me3_ac_combine.shape)
    # print(WT_me3_ac_combine.head())
    #
    # HS_me3 = pd.read_csv(src_dir / f'CDS002D13C.e500.juice', header=None, sep='\s+')
    # print(HS_me3.shape)
    # print(HS_me3.head())
    #
    #
    # HS_ac = pd.read_csv(src_dir / f'CDS005D11DC.e500.juice', header=None, sep='\s+')
    # print(HS_ac.shape)
    # print(HS_ac.head())
    #
    # HS_me3_ac_combine = pd.concat([HS_me3, HS_ac]).sort_values(by=[1, 5])
    # print(HS_me3_ac_combine.shape)
    # print(HS_me3_ac_combine.head())
    #
    # HS_me3_ac_combine.to_csv(src_dir / f'HS_me3_ac_combined.e500.juice', header=None, index=None, sep=' ')
    # HS_me3_ac_combine = pd.read_csv(src_dir / f'HS_me3_ac_combined.e500.juice', header=None, sep='\s+',
    #                                 dtype={1: str, 5: str})
    #
    # print(HS_me3_ac_combine.shape)
    # print(HS_me3_ac_combine.head())

