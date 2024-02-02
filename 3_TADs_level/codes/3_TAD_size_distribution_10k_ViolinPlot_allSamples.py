import pathlib as p
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from matplotlib.colors import LinearSegmentedColormap

# from scipy.stats import pearsonr, spearmanr, gaussian_kde

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    dest_dir = root.parent / 'results' / 'home_data_results' / '3_results_TAD_size_distribution_VioloinPlot'
    dest_dir.mkdir(parents=True, exist_ok=True)

    WT_num = '01'
    WT_color = '#4D7731'
    # WT_num = '04'
    # WT_color = '#F50CB8'

    src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    WT_TAD_raw = pd.read_csv(src_dir / f'CDS0{WT_num}D_domains.bedpe', header=None, sep='\t')
    WT_TAD = WT_TAD_raw.iloc[:, :3]
    WT_TAD.columns = ['chromosome', 'start', 'end']
    print(WT_TAD.shape)
    # print(WT_TAD.head())
    WT_TAD_addSize = WT_TAD.copy()
    WT_TAD_addSize['size'] = WT_TAD.apply(lambda x: x['end'] -x['start'], axis=1) / 1000
    # print(WT_TAD_addSize.shape)
    # print(WT_TAD_addSize.head())
    print(WT_TAD_addSize['size'].mean())
    print(WT_TAD_addSize['size'].median())
    WT_min = WT_TAD_addSize['size'].min()
    WT_max = WT_TAD_addSize['size'].max()
    df1 = WT_TAD_addSize['size'].reset_index()
    df1['Category'] = 'H3K27me3_WT'
    df1.drop(columns='index', inplace=True)

    HS_num = '02'
    HS_color = '#98A246'
    # HS_num = '05'
    # HS_color = '#743CA6'
    src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    HS_TAD_raw = pd.read_csv(src_dir / f'CDS0{HS_num}D_domains.bedpe', header=None, sep='\t')
    HS_TAD = HS_TAD_raw.iloc[:, :3]
    HS_TAD.columns = ['chromosome', 'start', 'end']
    print(HS_TAD.shape)
    # print(HS_TAD.head())
    HS_TAD_addSize = HS_TAD.copy()
    HS_TAD_addSize['size'] = HS_TAD.apply(lambda x: x['end'] -x['start'], axis=1) / 1000
    # print(HS_TAD_addSize.shape)
    # print(HS_TAD_addSize.head())
    print(HS_TAD_addSize['size'].mean())
    print(HS_TAD_addSize['size'].median())
    HS_min = HS_TAD_addSize['size'].min()
    HS_max = HS_TAD_addSize['size'].max()
    df2 = HS_TAD_addSize['size'].reset_index()
    df2['Category'] = 'H3K27me3_HS'
    df2.drop(columns='index', inplace=True)

    #
    ac_WT_num = '04'
    ac_WT_color = '#F50CB8'
    src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{ac_WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    ac_WT_TAD_raw = pd.read_csv(src_dir / f'CDS0{ac_WT_num}D_domains.bedpe', header=None, sep='\t')
    ac_WT_TAD = ac_WT_TAD_raw.iloc[:, :3]
    ac_WT_TAD.columns = ['chromosome', 'start', 'end']
    print(ac_WT_TAD.shape)
    # print(ac_WT_TAD.head())
    ac_WT_TAD_addSize = ac_WT_TAD.copy()
    ac_WT_TAD_addSize['size'] = ac_WT_TAD.apply(lambda x: x['end'] -x['start'], axis=1) / 1000
    # print(ac_WT_TAD_addSize.shape)
    # print(ac_WT_TAD_addSize.head())
    print(ac_WT_TAD_addSize['size'].mean())
    print(ac_WT_TAD_addSize['size'].median())
    ac_WT_min = ac_WT_TAD_addSize['size'].min()
    ac_WT_max = ac_WT_TAD_addSize['size'].max()
    df3 = ac_WT_TAD_addSize['size'].reset_index()
    df3['Category'] = 'H3K27ac_WT'
    df3.drop(columns='index', inplace=True)


    ac_HS_num = '05'
    ac_HS_color = '#743CA6'
    src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{ac_HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    ac_HS_TAD_raw = pd.read_csv(src_dir / f'CDS0{ac_HS_num}D_domains.bedpe', header=None, sep='\t')
    ac_HS_TAD = ac_HS_TAD_raw.iloc[:, :3]
    ac_HS_TAD.columns = ['chromosome', 'start', 'end']
    print(ac_HS_TAD.shape)
    # print(ac_HS_TAD.head())
    ac_HS_TAD_addSize = ac_HS_TAD.copy()
    ac_HS_TAD_addSize['size'] = ac_HS_TAD.apply(lambda x: x['end'] -x['start'], axis=1) / 1000
    # print(ac_HS_TAD_addSize.shape)
    # print(ac_HS_TAD_addSize.head())
    print(ac_HS_TAD_addSize['size'].mean())
    print(ac_HS_TAD_addSize['size'].median())
    ac_HS_min = ac_HS_TAD_addSize['size'].min()
    ac_HS_max = ac_HS_TAD_addSize['size'].max()
    df4 = ac_HS_TAD_addSize['size'].reset_index()
    df4['Category'] = 'H3K27ac_HS'
    df4.drop(columns='index', inplace=True)

    combined_df = pd.concat([df1, df2, df3, df4])
    # print(combined_df)

    colors = {'H3K27me3_WT': '#4D7731', 'H3K27me3_HS': '#98A246', 'H3K27ac_WT': '#F50CB8', 'H3K27ac_HS': '#743CA6'}

    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(1)
    sns.violinplot(x='Category', y='size', data=combined_df, palette=colors,linewidth=1)
    ax.set_ylabel("TAD size (Kb)")
    ax.set_xlabel("")

    # plt.show()
    fig.savefig(dest_dir / f'{WT_num}_{HS_num}_{ac_WT_num}_{ac_HS_num}_ViolinPlot_TAD_size_allSamples.png', dpi=300, bbox_inches='tight')