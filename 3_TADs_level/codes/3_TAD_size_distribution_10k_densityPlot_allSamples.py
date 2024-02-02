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
    dest_dir = root.parent / 'results' / 'home_data_results' / '3_results_TAD_size_distribution'
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
    print(WT_TAD.head())
    WT_TAD_addSize = WT_TAD.copy()
    WT_TAD_addSize['size'] = WT_TAD.apply(lambda x: x['end'] -x['start'], axis=1) / 1000
    print(WT_TAD_addSize.shape)
    print(WT_TAD_addSize.head())
    print(WT_TAD_addSize['size'].mean())
    print(WT_TAD_addSize['size'].median())
    WT_min = WT_TAD_addSize['size'].min()
    WT_max = WT_TAD_addSize['size'].max()

    HS_num = '02'
    HS_color = '#98A246'
    # HS_num = '05'
    # HS_color = '#743CA6'
    src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    HS_TAD_raw = pd.read_csv(src_dir / f'CDS0{HS_num}D_domains.bedpe', header=None, sep='\t')
    HS_TAD = HS_TAD_raw.iloc[:, :3]
    HS_TAD.columns = ['chromosome', 'start', 'end']
    print(HS_TAD.shape)
    print(HS_TAD.head())
    HS_TAD_addSize = HS_TAD.copy()
    HS_TAD_addSize['size'] = HS_TAD.apply(lambda x: x['end'] -x['start'], axis=1) / 1000
    print(HS_TAD_addSize.shape)
    print(HS_TAD_addSize.head())
    print(HS_TAD_addSize['size'].mean())
    print(HS_TAD_addSize['size'].median())
    HS_min = HS_TAD_addSize['size'].min()
    HS_max = HS_TAD_addSize['size'].max()
    
    
    ac_WT_num = '04'
    ac_WT_color = '#F50CB8'
    src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{ac_WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    ac_WT_TAD_raw = pd.read_csv(src_dir / f'CDS0{ac_WT_num}D_domains.bedpe', header=None, sep='\t')
    ac_WT_TAD = ac_WT_TAD_raw.iloc[:, :3]
    ac_WT_TAD.columns = ['chromosome', 'start', 'end']
    print(ac_WT_TAD.shape)
    print(ac_WT_TAD.head())
    ac_WT_TAD_addSize = ac_WT_TAD.copy()
    ac_WT_TAD_addSize['size'] = ac_WT_TAD.apply(lambda x: x['end'] -x['start'], axis=1) / 1000
    print(ac_WT_TAD_addSize.shape)
    print(ac_WT_TAD_addSize.head())
    print(ac_WT_TAD_addSize['size'].mean())
    print(ac_WT_TAD_addSize['size'].median())
    ac_WT_min = ac_WT_TAD_addSize['size'].min()
    ac_WT_max = ac_WT_TAD_addSize['size'].max()


    ac_HS_num = '05'
    ac_HS_color = '#743CA6'
    src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{ac_HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    ac_HS_TAD_raw = pd.read_csv(src_dir / f'CDS0{ac_HS_num}D_domains.bedpe', header=None, sep='\t')
    ac_HS_TAD = ac_HS_TAD_raw.iloc[:, :3]
    ac_HS_TAD.columns = ['chromosome', 'start', 'end']
    print(ac_HS_TAD.shape)
    print(ac_HS_TAD.head())
    ac_HS_TAD_addSize = ac_HS_TAD.copy()
    ac_HS_TAD_addSize['size'] = ac_HS_TAD.apply(lambda x: x['end'] -x['start'], axis=1) / 1000
    print(ac_HS_TAD_addSize.shape)
    print(ac_HS_TAD_addSize.head())
    print(ac_HS_TAD_addSize['size'].mean())
    print(ac_HS_TAD_addSize['size'].median())
    ac_HS_min = ac_HS_TAD_addSize['size'].min()
    ac_HS_max = ac_HS_TAD_addSize['size'].max()

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(1)
    sns.kdeplot(WT_TAD_addSize['size'], ax=ax, lw=1, color=WT_color)
    sns.kdeplot(HS_TAD_addSize['size'], ax=ax, lw=1, color=HS_color)
    sns.kdeplot(ac_WT_TAD_addSize['size'], ax=ax, lw=1, color=ac_WT_color)
    sns.kdeplot(ac_HS_TAD_addSize['size'], ax=ax, lw=1, color=ac_HS_color)
    # ax.axvline(x=100, color='gray', linestyle='--')
    # ax.axvline(x=140, color='gray', linestyle='--')
    ax.set_ylabel('Scaled density')
    ax.set_xticks([200, 600, 1000])
    ax.set_xlabel('TAD size (Kb)')
    plt.show()
    # fig.savefig(dest_dir / f'{WT_num}_{HS_num}_{ac_WT_num}_{ac_HS_num}_densityPlot_TAD_size_allSamples.png', dpi=300, bbox_inches='tight')