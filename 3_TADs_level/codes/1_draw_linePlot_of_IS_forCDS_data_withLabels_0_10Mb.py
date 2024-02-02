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

# from scipy.stats import pearsonr, spearmanr, gaussian_kde

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    # chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX']
    # for chro in chro_list:
    #     print(chro)
    chro = '2L'

    WT_num = '04'
    HS_num = '05'
    # WT_color = '#4D7731'
    # HS_color = '#98A246'

    WT_color = '#F50CB8'
    HS_color = '#743CA6'
    name = 'H3K27ac'
    x=0
    y=1000
    # x = 30
    # y = 320

    WT_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{WT_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    HS_src_dir = root.parent / 'results' / 'home_data_results' / f'CDS0{HS_num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh0.01_delta0.01_dfr'
    dest_dir = root.parent / 'results' / 'home_data_results' / f'1_results_linePlot_of_IS_forCDS_data'
    dest_dir.mkdir(parents=True, exist_ok=True)

    WT_IS = pd.read_csv(WT_src_dir / f'CDS0{WT_num}D_score.bedgraph', header=None, sep='\t')
    WT_IS.columns = ['chromosome', 'start', 'end', 'score']
    WT_IS_addMidpoints = WT_IS.copy()
    WT_IS_addMidpoints['mid'] = WT_IS.apply(lambda x: 0.5*(x['start'] + x['end']), axis=1)
    WT_IS_addMidpoints['mid'] = WT_IS_addMidpoints['mid'] / 10000
    WT_IS_addMidpoints_fiterChom = WT_IS_addMidpoints[WT_IS_addMidpoints['chromosome'] == chro].reset_index(drop=True).iloc[x:y]
    # print(WT_IS_addMidpoints_fiterChom)
    # print(WT_IS_addMidpoints_fiterChom.shape)
    # print(WT_IS_addMidpoints_fiterChom.head(10))
    # print(WT_IS_addMidpoints_fiterChom[['score', 'mid']])

    HS_IS = pd.read_csv(HS_src_dir / f'CDS0{HS_num}D_score.bedgraph', header=None, sep='\t')
    HS_IS.columns = ['chromosome', 'start', 'end', 'score']
    HS_IS_addMidpoints = HS_IS.copy()
    HS_IS_addMidpoints['mid'] = HS_IS.apply(lambda x: 0.5 * (x['start'] + x['end']), axis=1)
    HS_IS_addMidpoints['mid'] = HS_IS_addMidpoints['mid'] / 10000
    HS_IS_addMidpoints_fiterChom = HS_IS_addMidpoints[HS_IS_addMidpoints['chromosome'] == chro].reset_index(drop=True).iloc[x:y]


    fig, ax = plt.subplots(figsize=(8, 0.8), dpi=300)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(0.5)

    plt.plot(WT_IS_addMidpoints_fiterChom['mid'], WT_IS_addMidpoints_fiterChom['score'], color=WT_color, linewidth=0.5)
    plt.plot(HS_IS_addMidpoints_fiterChom['mid'], HS_IS_addMidpoints_fiterChom['score'], color=HS_color, linewidth=0.5)
    # plt.axhline(y=0, color='black', linestyle='--', linewidth=0.5)


    #
    ticks = [i for i in range(x, y, 100)]
    print(ticks)
    # labels = ['0', '100K', '200K', '300K', '400K', '500K', '600K', '700K', '800K', '900K', '1M']
    labels = ['0', '1M', '2M', '3M', '4M', '5M', '6M', '7M', '8M', '9M', '10M']
    # ax.set_xticks(ticks)
    ax.set_xticks([])
    ax.set_xlim(x, y)
    # ax.set_xticklabels('')
    # ax.set_xticklabels(labels)
    y_ticks = [-1, 0 , 1]
    ax.set_yticks([])
    # ax.set_yticks(y_ticks)
    # ax.set_yticklabels(y_ticks)
    # ax.set_yticklabels('')
    #
    # plt.xlabel('')
    # plt.ylabel('')
    # plt.title('')
    # plt.show()

    plt.savefig(dest_dir / f'linePlot_with_correlation_of_IS_vector_{name}_{chro}_WT_HS_10k_0_10Mb.png', dpi=300, bbox_inches='tight')

    # groups = [['01', '02', '#4D7731', '#98A246', 'H3K27me3'], ['04', '05', '#F50CB8', '#743CA6', 'H3K27ac'], ['07', '08', '#1087F4', '#99BDCB', 'RNAP2']]
    # for group in groups:
    #     WT_num = group[0]
    #     HS_num = group[1]
    #     WT_color = group[2]
    #     HS_color = group[3]
    #     name = group[4]