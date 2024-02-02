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
    dest_dir = root.parent / 'results' / 'home_data_results' / '3_results_TAD_size_distribution_withoutLabel'
    dest_dir.mkdir(parents=True, exist_ok=True)

    # WT_num = '01'
    # WT_color = '#4D7731'
    WT_num = '04'
    WT_color = '#F50CB8'

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
    print(WT_TAD_addSize['size'].max())
    WT_min = WT_TAD_addSize['size'].min()
    WT_max = WT_TAD_addSize['size'].max()

    # HS_num = '02'
    # HS_color = '#98A246'
    HS_num = '05'
    HS_color = '#743CA6'
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
    print(HS_TAD_addSize['size'].max())
    HS_min = HS_TAD_addSize['size'].min()
    HS_max = HS_TAD_addSize['size'].max()

    bin_edges = np.arange(0, 1000, 100)
    hist1, _ = np.histogram(WT_TAD_addSize['size'], bins=bin_edges)
    hist2, _ = np.histogram(HS_TAD_addSize['size'], bins=bin_edges)

    # Set width for the bars
    width = 40  # adjust width as needed

    bar_positions1 = bin_edges[:-1] - width / 2  # for series1
    bar_positions2 = bin_edges[:-1] + width / 2

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.bar(bar_positions1, hist1, width=width, align='center', color=WT_color)
    plt.bar(bar_positions2, hist2, width=width, align='center', color=HS_color)

    # plt.xticks(bin_edges[:-1], [f"{edge}-{bin_edges[i + 1]}" for i, edge in enumerate(bin_edges[:-1])], rotation=45)
    ax.set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800])
    ax.set_xticklabels([])
    ax.set_yticks([100, 200, 300, 400])
    ax.set_yticklabels([])
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    plt.show()
    # fig.savefig(dest_dir / f'{WT_num}_{HS_num}_histogram_TAD_size_withoutLabel.png', dpi=300, bbox_inches='tight')