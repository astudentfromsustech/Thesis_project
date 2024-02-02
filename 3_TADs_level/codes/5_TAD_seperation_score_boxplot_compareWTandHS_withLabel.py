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

    # groups = [['01', '02', 'H3K27me3'], ['04', '05', 'H3K27ac'], ['07', '08', 'RNAP2']]
    groups = [['04', '05', 'H3K27ac']]
    for group in groups:
        WT_num = group[0]
        HS_num = group[1]
        name = group[2]
        print(WT_num)
        print(HS_num)
        print(name)

        WT_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
        HS_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
        dest_dir = root.parent / 'results' / 'home_data_results' / f'5_TAD_seperation_score_boxplot_withDots_jitter_withLabel'
        dest_dir.mkdir(parents=True, exist_ok=True)

        chro_list = ['2L', '2R', '3L', '3R', 'X']
        # chro_list = ['2L']
        for chro in chro_list:
            print(chro)

            WT_IS = pd.read_csv(WT_src_dir / f'CDS0{WT_num}D_boundaries.bed', header=None, sep='\t')
            WT_IS.columns = ['chromosome', 'start', 'end', 'ID', 'score', 'dot']
            # print(WT_IS.shape)
            # print(WT_IS.head())

            WT_IS_score = WT_IS[WT_IS['chromosome'] == chro]['score']
            print(WT_IS_score.mean())

            HS_IS = pd.read_csv(HS_src_dir / f'CDS0{HS_num}D_boundaries.bed', header=None, sep='\t')
            HS_IS.columns = ['chromosome', 'start', 'end', 'ID', 'score', 'dot']
            # print(HS_IS.shape)
            # print(HS_IS.head())

            HS_IS_score = HS_IS[HS_IS['chromosome'] == chro]['score']
            print(HS_IS_score.mean())

            # t_stat, p_val= ttest_ind(WT_IS_score,  HS_IS_score)
            m_stat, p_val = mannwhitneyu(WT_IS_score,  HS_IS_score, alternative='less')

            # print("T-statistic:", t_stat)
            print("P-value:", p_val)

            df = pd.DataFrame(
                {'WT': WT_IS_score, 'HS': HS_IS_score})
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

            ax.set_title(f"{name}_{chro}_WT_HS\nmannwhitneyu_pvalue_{p_val:.2f}", fontsize=8)
            plt.show()
            # plt.savefig(dest_dir / f'{name}_{chro}_WT_HS_mannwhitneyu_pvalue_{p_val:.2f}_10k_withLabel.png', dpi=300, bbox_inches='tight')

