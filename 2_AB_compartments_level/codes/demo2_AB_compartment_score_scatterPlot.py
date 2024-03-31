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

from scipy.stats import pearsonr, spearmanr, gaussian_kde, ttest_ind

import sys
print(sys.float_info.min)


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '1_ABcompartments_classification'
    dest_dir = root.parent / 'results' / 'demo2_AB_compartment_score_scatterPlot'
    dest_dir.mkdir(parents=True, exist_ok=True)

    resolution = 10000

    WT = pd.read_csv(src_dir / f'WT_compartment_score_resolution{resolution}.bedGraph', header=None, sep='\t')
    WT.columns = ['chromosome', 'start', 'end', 'compartment_score', 'compartment']
    # print(WT.shape)
    # print(WT.head())
    WT.loc[WT['compartment'] == 'A', 'compartment_score'] = WT.loc[WT['compartment'] == 'A', 'compartment_score'] * (-1)
    WT.loc[WT['compartment'] == 'B', 'compartment_score'] = WT.loc[WT['compartment'] == 'B', 'compartment_score'] * (-1)

    # Filter WT by compartment after modification
    WT_A = WT[WT['compartment'] == 'A']
    WT_B = WT[WT['compartment'] == 'B']

    print(WT_A.shape)
    print(WT_A.head())
    print(WT_B.shape)
    print(WT_B.head())

    # Reading HS compartment score data
    HS = pd.read_csv(src_dir / f'HS_compartment_score_resolution{resolution}.bedGraph', header=None, sep='\t')
    HS.columns = ['chromosome', 'start', 'end', 'compartment_score', 'compartment']
    # print(HS.shape)
    # print(HS.head())
    # Modifying HS compartment scores directly in the DataFrame
    HS.loc[HS['compartment'] == 'A', 'compartment_score'] = HS.loc[HS['compartment'] == 'A', 'compartment_score'] * (-1)
    HS.loc[HS['compartment'] == 'B', 'compartment_score'] = HS.loc[HS['compartment'] == 'B', 'compartment_score'] * (-1)
    # Filter HS by compartment after modification
    HS_A = HS[HS['compartment'] == 'A']
    HS_B = HS[HS['compartment'] == 'B']
    print(HS_A.shape)
    print(HS_A.head())
    print(HS_B.shape)
    print(HS_B.head())

    t_stat_pos, p_value_pos = ttest_ind(WT_A['compartment_score'], HS_A['compartment_score'])

    print(f"T-statistic: {t_stat_pos}")
    print(f"P-value: {p_value_pos}")

    t_stat_neg, p_value_neg = ttest_ind(WT_B['compartment_score'], HS_B['compartment_score'])

    print(f"T-statistic: {t_stat_neg}")
    print(f"P-value: {p_value_neg}")


    df = pd.DataFrame({'WT_A': WT_A['compartment_score'], 'HS_A': HS_A['compartment_score'], 'WT_B': WT_B['compartment_score'], 'HS_B': HS_B['compartment_score']})
    # print(df)
    fig, ax = plt.subplots(figsize=(5, 5), dpi=200)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(1)
    # # #
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    colors_pal = ['#ed8467', '#f1cab6', '#7b9ef8', '#c0d3f5']
    flierprops = dict(marker='o', markerfacecolor='none', markersize=4, markeredgecolor='black', linestyle='none',
                      markeredgewidth=0.6)
    sns.boxplot(data=df, palette=colors_pal, flierprops=flierprops)
    sns.stripplot(data=df, jitter=True, marker='o', alpha=0.7, color='black', s=1)
    # ax.set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
    # ax.set_yticklabels([])
    # ax.set_xticklabels([])

    # ax.set_title(f"{name}_{chro}_WT_HS\nA_pvalue_{p_value_pos:.2f}_B_pvalue_{p_value_neg:.2f}", fontsize=8)
    # plt.savefig(dest_dir / f'{name}_{chro}_WT_HS_A_pvalue_{p_value_pos:.2f}_B_pvalue_{p_value_neg:.2f}_50k_KR_withoutLabel.png', dpi=300, bbox_inches='tight')
    plt.show()
