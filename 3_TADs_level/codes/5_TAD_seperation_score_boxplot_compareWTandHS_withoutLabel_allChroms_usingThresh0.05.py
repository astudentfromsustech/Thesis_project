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

from scipy.stats import pearsonr, spearmanr, gaussian_kde, ttest_ind, mannwhitneyu

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent

    groups = [['01', '02', 'H3K27me3'], ['04', '05', 'H3K27ac'], ['07', '08', 'RNAP2']]
    # groups = [['01', '02', 'H3K27me3']]
    for group in groups:
        WT_num = group[0]
        HS_num = group[1]
        name = group[2]
        print(WT_num)
        print(HS_num)
        print(name)

        WT_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.05_delta0.01_dfr'
        HS_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.05_delta0.01_dfr'
        dest_dir = root.parent / 'results' / 'home_data_results' / f'5_TAD_seperation_score_boxplot_withDots_jitter_withoutLabel_allChroms_usingThresh0.05'
        dest_dir.mkdir(parents=True, exist_ok=True)

        WT_IS_score_allChroms = pd.Series(dtype='float64')
        HS_IS_score_allChroms = pd.Series(dtype='float64')
        chro_list = ['2L', '2R', '3L', '3R', 'X']
        # chro_list = ['2L']
        for chro in chro_list:
            print(chro)

            WT_IS = pd.read_csv(WT_src_dir / f'CDS0{WT_num}D_boundaries.bed', header=None, sep='\t')
            WT_IS.columns = ['chromosome', 'start', 'end', 'ID', 'score', 'dot']
            # print(WT_IS.shape)
            # print(WT_IS.head())

            WT_IS_score = WT_IS[WT_IS['chromosome'] == chro]['score']


            HS_IS = pd.read_csv(HS_src_dir / f'CDS0{HS_num}D_boundaries.bed', header=None, sep='\t')
            HS_IS.columns = ['chromosome', 'start', 'end', 'ID', 'score', 'dot']
            # print(HS_IS.shape)
            # print(HS_IS.head())

            HS_IS_score = HS_IS[HS_IS['chromosome'] == chro]['score']


            # t_stat, p_val= ttest_ind(WT_IS_score,  HS_IS_score)

            WT_IS_score_allChroms = pd.concat([WT_IS_score_allChroms, WT_IS_score]).reset_index(drop=True)
            HS_IS_score_allChroms = pd.concat([HS_IS_score_allChroms, HS_IS_score]).reset_index(drop=True)

        df = pd.DataFrame(
            {'WT': WT_IS_score_allChroms, 'HS': HS_IS_score_allChroms})

        m_stat, p_val = mannwhitneyu(WT_IS_score, HS_IS_score, alternative='less')
        print(WT_IS_score_allChroms.mean())
        print(HS_IS_score_allChroms.mean())
        # print("T-statistic:", t_stat)
        print("P-value:", p_val)
        fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(1)
        # # #

        colors_pal = ['#90EE90', '#f1cab6']
        flierprops = dict(marker='o', markerfacecolor='none', markersize=4, markeredgecolor='black',
                          linestyle='none', markeredgewidth=0.6)
        sns.boxplot(data=df, palette=colors_pal, flierprops=flierprops)
        sns.stripplot(data=df, jitter=True, marker='o', alpha=0.7, color='black', s=1)
        ax.set_xticklabels([])
        ax.set_yticks([-1.5, -1, -0.5, 0, 0.5])
        ax.set_yticklabels([])
        # ax.set_title(f"{name}_allChroms_WT_HS\nmannwhitneyu_pvalue_{p_val:.2f}", fontsize=8)
        # plt.show()
        plt.savefig(dest_dir / f'{name}_allChroms_WT_HS_mannwhitneyu_pvalue_{p_val}_10k_withoutLabel_allChroms.png', dpi=300, bbox_inches='tight')

