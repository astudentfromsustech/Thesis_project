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
    dest_dir = root.parent / 'results' / 'PNAS_data_V2_results' / '3_results_TAD_size_distribution_forPNAS_S2'
    dest_dir.mkdir(parents=True, exist_ok=True)


    WT_color = 'black'

    src_dir = root.parent / 'results' / 'PNAS_data_V2_results' / f'WT_combined_10000_targetChroms_normalize0_1_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    WT_TAD_raw = pd.read_csv(src_dir / f'WT_domains.bedpe', header=None, sep='\t')
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


    HS_color = 'red'
    src_dir = root.parent / 'results' / 'PNAS_data_V2_results' / f'HS_combined_10000_targetChroms_normalize0_1_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    HS_TAD_raw = pd.read_csv(src_dir / f'HS_domains.bedpe', header=None, sep='\t')
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
    ax.axvline(x=100, color='gray', linestyle='--')
    ax.set_ylabel('Scaled density')
    ax.set_xticks([200, 600, 1000])
    ax.set_xlabel('TAD size (Kb)')
    plt.show()
    # fig.savefig(dest_dir / f'WT_HS_scaled_density_TAD_size_PNAS_S2.png', dpi=300, bbox_inches='tight')