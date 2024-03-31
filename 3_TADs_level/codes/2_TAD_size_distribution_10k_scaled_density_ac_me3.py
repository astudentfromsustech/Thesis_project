# import cooler
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
    dest_dir = root.parent / 'results' / 'home_data_results' / '2_TAD_size_distribution_10k_histogram_for_ac_me3'
    dest_dir.mkdir(parents=True, exist_ok=True)

    # WT_num = '01'
    # WT_color = '#4D7731'
    WT_num = 'WT'
    WT_color = '#4D7731'
    # WT_num = '04'
    # WT_color = '#F50CB8'

    src_dir = root.parent / 'results' / f'ac_me3_downsample_{WT_num}_me3_ac_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.1_delta0.01_dfr'
    WT_TAD_raw = pd.read_csv(src_dir / f'ac_me3_downsample_{WT_num}_me3_ac_combined_res10000_domains.bedpe', header=None, sep='\t')
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
    WT_kde = gaussian_kde(WT_TAD_addSize['size'])
    WT_x_range = np.linspace(WT_min, WT_max, 1000)
    WT_y_vals = WT_kde(WT_x_range)
    WT_y_vals = WT_y_vals / WT_y_vals.max()

    # HS_num = '02'
    # HS_color = '#98A246'
    HS_num = 'HS'
    HS_color = '#98A246'
    # HS_num = '05'
    # HS_color = '#743CA6'
    src_dir = root.parent / 'results' / f'ac_me3_downsample_{HS_num}_me3_ac_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.1_delta0.01_dfr'
    HS_TAD_raw = pd.read_csv(src_dir / f'ac_me3_downsample_{HS_num}_me3_ac_combined_res10000_domains.bedpe',
                             header=None, sep='\t')
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
    HS_kde = gaussian_kde(HS_TAD_addSize['size'])
    HS_x_range = np.linspace(HS_min, HS_max, 1000)
    HS_y_vals = HS_kde(HS_x_range)
    HS_y_vals = HS_y_vals / HS_y_vals.max()


    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(1)
    plt.plot(WT_x_range, WT_y_vals, lw=1, color=WT_color)
    plt.plot(HS_x_range, HS_y_vals, lw=1, color=HS_color)
    # ax.axvline(x=100, color='gray', linestyle='--')

    ax.set_xticks([200, 600, 1000])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.show()
    # fig.savefig(dest_dir / f'{WT_num}_{HS_num}_scaled_density_TAD_size_withoutLabel.png', dpi=300, bbox_inches='tight')