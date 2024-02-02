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
    # chro_list = ['2L', '2R', '3L', '3R', 'X']
    chro_list = ['2L']
    for chro in chro_list:
        print(chro)
        # WT_num = '01'
        # HS_num = '02'
        WT_num = '04'
        HS_num = '05'

        WT_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
        HS_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
        dest_dir = root.parent / 'results' / 'home_data_results' / f'1_results_ScatterPlot_withColorDensity_of_IS_forCDS_data'
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

        x = result['score_x']
        y = result['score_y']
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)

        # Sort the points by density, so the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

        fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            # spine.set_color('none')
            spine.set_linewidth(1)

        sc = ax.scatter(x, y, c=z, s=1)
        cbar=plt.colorbar(sc, ax=ax, label='Density')

        ax.set_xticks([-1,0,1,2,3])
        ax.set_yticks([-1, 0, 1,2,3])
        cbar.set_ticks([])
        cbar.set_label('')
        plt.show()
        # plt.savefig(dest_dir / f'ScatterPlot_with_correlation_of_IS_vector_CDS0{WT_num}D_CDS0{HS_num}D_{chro}_pearson_{round(pearson_corr, 2)}_spearman_{round(spearman_corr, 2)}.png', dpi=300,
        #             bbox_inches='tight')

