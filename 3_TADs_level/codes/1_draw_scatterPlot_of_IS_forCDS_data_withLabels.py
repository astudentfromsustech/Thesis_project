# import cooler
import pathlib as p
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from matplotlib.colors import LinearSegmentedColormap

from scipy.stats import pearsonr, spearmanr, gaussian_kde

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    # chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX']
    # for chro in chro_list:
    #     print(chro)
    chro = '2L'

    WT_num = '01'
    HS_num = '02'

    name = 'H3K27me3'


    WT_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    HS_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    dest_dir = root.parent / 'results' / 'home_data_results' / f'1_results_ScatterPlot_of_IS_forCDS_data'
    dest_dir.mkdir(parents=True, exist_ok=True)

    WT_IS = pd.read_csv(WT_src_dir / f'CDS0{WT_num}D_score.bedgraph', header=None, sep='\t')
    WT_IS.columns = ['chromosome', 'start', 'end', 'score']
    WT_IS_addMidpoints = WT_IS.copy()
    WT_IS_addMidpoints['mid'] = WT_IS.apply(lambda x: 0.5*(x['start'] + x['end']), axis=1)
    WT_IS_addMidpoints['mid'] = WT_IS_addMidpoints['mid'] / 10000
    WT_IS_addMidpoints_fiterChom = WT_IS_addMidpoints[WT_IS_addMidpoints['chromosome'] == chro].reset_index(drop=True)
    print(WT_IS_addMidpoints_fiterChom.shape)
    # print(WT_IS_addMidpoints_fiterChom)
    # print(WT_IS_addMidpoints_fiterChom.shape)
    # print(WT_IS_addMidpoints_fiterChom.head(10))
    # print(WT_IS_addMidpoints_fiterChom[['score', 'mid']])

    HS_IS = pd.read_csv(HS_src_dir / f'CDS0{HS_num}D_score.bedgraph', header=None, sep='\t')
    HS_IS.columns = ['chromosome', 'start', 'end', 'score']
    HS_IS_addMidpoints = HS_IS.copy()
    HS_IS_addMidpoints['mid'] = HS_IS.apply(lambda x: 0.5 * (x['start'] + x['end']), axis=1)
    HS_IS_addMidpoints['mid'] = HS_IS_addMidpoints['mid'] / 10000
    HS_IS_addMidpoints_fiterChom = HS_IS_addMidpoints[HS_IS_addMidpoints['chromosome'] == chro].reset_index(drop=True)
    print(HS_IS_addMidpoints_fiterChom.shape)

    result = WT_IS_addMidpoints_fiterChom.merge(HS_IS_addMidpoints_fiterChom, on=['chromosome','start', 'end', 'mid'], how='inner')
    print(result.shape)
    print(result.tail())

    pearson_corr, pearson_pvalue = pearsonr(result['score_x'], result['score_y'])

    # Calculate Spearman's correlation using scipy
    spearman_corr, spearman_pvalue = spearmanr(result['score_x'], result['score_y'])
    #
    print("Pearson's correlation:", pearson_corr)
    print("Pearson's p-value:", pearson_pvalue)

    print("\nSpearman's correlation:", spearman_corr)
    print("Spearman's p-value:", spearman_pvalue)

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(0.5)
    plt.scatter(result['score_x'], result['score_y'], s=1, color='black', alpha=0.7)
    plt.show()
    #
    # plt.plot(WT_IS_addMidpoints_fiterChom['mid'], WT_IS_addMidpoints_fiterChom['score'], color=WT_color, linewidth=0.5)
    # plt.plot(HS_IS_addMidpoints_fiterChom['mid'], HS_IS_addMidpoints_fiterChom['score'], color=HS_color, linewidth=0.5)
    # plt.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    # # plt.axvline(x=1, color='black', linestyle='-', linewidth=1)  # red solid line at x=1
    # # plt.axvline(x=1000, color='black', linestyle='-', linewidth=1)
    #
    # #
    # ticks = [i for i in range(x, y, 100)]
    # print(ticks)
    # # labels = ['0', '100K', '200K', '300K', '400K', '500K', '600K', '700K', '800K', '900K', '1M']
    # labels = ['0', '1M', '2M', '3M', '4M', '5M', '6M', '7M', '8M', '9M', '10M']
    # # ax.set_xticks(ticks)
    # ax.set_xticks([])
    # ax.set_xlim(x, y)
    # # ax.set_xticklabels('')
    # # ax.set_xticklabels(labels)
    # y_ticks = [-1, 0 , 1]
    # ax.set_yticks([])
    # # ax.set_yticks(y_ticks)
    # # ax.set_yticklabels(y_ticks)
    # # ax.set_yticklabels('')
    # #
    # # plt.xlabel('')
    # # plt.ylabel('')
    # # plt.title('')
    # # plt.show()
    #
    # # plt.savefig(dest_dir / f'ScatterPlot_with_correlation_of_IS_vector_{name}_{chro}_WT_HS_10k_300k-3Mb.png', dpi=300, bbox_inches='tight')
