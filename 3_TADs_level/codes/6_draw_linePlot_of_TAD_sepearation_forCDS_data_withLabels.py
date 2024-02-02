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

from statsmodels.nonparametric.smoothers_lowess import lowess
# from scipy.stats import pearsonr, spearmanr, gaussian_kde

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def moving_average(data, window_size):
    """Compute a simple moving average."""
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    #
    # chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX']
    # for chro in chro_list:
    #     print(chro)

    # chro = 'chr2L'
    # groups = [['01', '02', 'H3K27me3'], ['04', '05', 'H3K27ac'], ['07', '08', 'RNAP2']]
    groups = [['01', '02', 'H3K27me3']]
    for group in groups:
        WT_num = group[0]
        HS_num = group[1]
        name = group[2]
        print(WT_num)
        print(HS_num)
        print(name)

        WT_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
        HS_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
        dest_dir = root.parent / 'results' / 'home_data_results' / f'6_results_linePlot_of_TAD_seperation_forCDS_data_withLabel'
        dest_dir.mkdir(parents=True, exist_ok=True)

        # chro_list = ['2L', '2R', '3L', '3R', 'X']
        chro_list = ['2L']
        for chro in chro_list:
            print(chro)

        WT_IS = pd.read_csv(WT_src_dir / f'CDS0{WT_num}D_boundaries.bed', header=None, sep='\t')
        WT_IS.columns = ['chromosome', 'start', 'end', 'ID', 'score', 'dot']
        WT_IS_addMid = WT_IS.copy()
        WT_IS_addMid['mid'] = WT_IS.apply(lambda x: 0.5*(x['end'] + x['start']), axis=1)
        WT_IS_addMid['mid'] = WT_IS_addMid['mid'] / 10000
        print(WT_IS_addMid.shape)
        print(WT_IS_addMid.head())
        # print(WT_IS_addMid[['mid', 'score']])

        HS_IS = pd.read_csv(HS_src_dir / f'CDS0{HS_num}D_boundaries.bed', header=None, sep='\t')
        HS_IS.columns = ['chromosome', 'start', 'end', 'ID', 'score', 'dot']
        HS_IS_addMid = HS_IS.copy()
        HS_IS_addMid['mid'] = HS_IS.apply(lambda x: 0.5 * (x['end'] + x['start']), axis=1)
        HS_IS_addMid['mid'] = HS_IS_addMid['mid'] / 10000
        print(HS_IS_addMid.shape)
        print(HS_IS_addMid.head())
        # print(HS_IS_addMid[['mid', 'score']])

        fig, ax = plt.subplots(figsize=(8, 0.8), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            # spine.set_color('none')
            spine.set_linewidth(0.5)
    #
        window_size = 3

        WT_x_vals = WT_IS_addMid['mid'][:60]
        WT_y_vals = WT_IS_addMid['score'][:60]
        WT_y_vals_smoothed = moving_average(WT_y_vals, window_size)
        # Adjust x values since the smoothed data will be shorter by window_size-1
        WT_x_vals_adjusted = WT_x_vals[(window_size - 1):]

        HS_x_vals = HS_IS_addMid['mid'][:50]
        HS_y_vals = HS_IS_addMid['score'][:50]
        HS_y_vals_smoothed = moving_average(HS_y_vals, window_size)
        # Adjust x values since the smoothed data will be shorter by window_size-1
        HS_x_vals_adjusted = HS_x_vals[(window_size - 1):]

        # Plotting
        plt.plot(WT_x_vals, WT_y_vals, color='k', linewidth=0.5, label='Original Data')
        plt.plot(HS_x_vals, HS_y_vals, color='gray', linewidth=0.5, label='Original Data')
        # plt.plot(WT_x_vals_adjusted, WT_y_vals_smoothed, color='r', linewidth=0.5, linestyle='--', label='Smoothed Curve')
        # plt.plot(HS_x_vals_adjusted, HS_y_vals_smoothed, color='g', linewidth=0.5, linestyle='--',
        #          label='Smoothed Curve')
        # # plt.legend()
        plt.show()


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
    # # plt.savefig(dest_dir / f'linePlot_with_correlation_of_IS_vector_{name}_{chro}_WT_HS_10k_0_10Mb.png', dpi=300, bbox_inches='tight')
