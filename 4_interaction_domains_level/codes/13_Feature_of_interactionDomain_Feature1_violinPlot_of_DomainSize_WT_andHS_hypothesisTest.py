import pathlib as p
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')

from scipy.stats import pearsonr, spearmanr, gaussian_kde, ttest_ind, mannwhitneyu, wilcoxon, ks_2samp

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)




if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '11_call_interaction_domains_STEP3_dropIntervalLen'
    dest_dir = root.parent / 'results' / '13_Feature_of_interactionDomain_Feature1_violinPlot_of_DomainSize_WT_andHS_hypothesisTest'
    dest_dir.mkdir(parents=True, exist_ok=True)

    num_list = [['01', '02', '#4D7731', '#98A246'], ['04', '05', '#F50CB8', '#743CA6'], ['07', '08', '#1087F4', '#99BDCB']]
    # num_list = [['01', '02', '#4D7731', '#98A246']]
    for num_pair in num_list:
        WT_num = num_pair[0]
        HS_num = num_pair[1]
        WT_color = num_pair[2]
        HS_color = num_pair[3]

        WT_domain_raw = pd.read_csv(src_dir / f'CDS0{WT_num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
        WT_domain_raw.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        # print(WT_domain_raw.shape)
        # print(WT_domain_raw.head())
        WT_domain = WT_domain_raw[['chromosome', 'start', 'end']]
        # print(WT_domain.shape)
        # print(WT_domain.head())
        WT_domain_addSize = WT_domain.copy()
        WT_domain_addSize['size'] = WT_domain.apply(lambda row: row['end'] - row['start'], axis=1) / 1000
        # print(WT_domain_addSize.shape)
        # print(WT_domain_addSize.head())

        HS_domain_raw = pd.read_csv(src_dir / f'CDS0{HS_num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
        HS_domain_raw.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        # print(HS_domain_raw.shape)
        # print(HS_domain_raw.head())
        HS_domain = HS_domain_raw[['chromosome', 'start', 'end']]
        # print(HS_domain.shape)
        # print(HS_domain.head())
        HS_domain_addSize = HS_domain.copy()
        HS_domain_addSize['size'] = HS_domain.apply(lambda row: row['end'] - row['start'], axis=1) / 1000
        # print(HS_domain_addSize.shape)
        # print(HS_domain_addSize.head())

        WT_domain_size = WT_domain_addSize['size']
        HS_domain_size = HS_domain_addSize['size']
        m_stat, p_val = mannwhitneyu(WT_domain_size, HS_domain_size, alternative='less')
        # m_stat, p_val = ks_2samp(WT_domain_size, HS_domain_size, alternative='less')
        print(p_val)

        HS_domain_size_filter = HS_domain_size[HS_domain_size <= 1000]
        # df = pd.DataFrame(
        #     {'WT': WT_domain_size, 'HS': HS_domain_size})
        df = pd.DataFrame({
            'Value': pd.concat([WT_domain_size, HS_domain_size_filter]),
            'Group': ['WT'] * len(WT_domain_size) + ['HS'] * len(HS_domain_size_filter)
        })

        color_palette = {"WT": WT_color, "HS": HS_color}
        fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            # spine.set_color('none')
            spine.set_linewidth(1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # colors_pal = [WT_color, HS_color]
        sns.violinplot(x='Group', y='Value', data=df, palette=color_palette)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([])
        # ax.set_yticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200])
        ax.set_yticks([0, 100, 200, 300, 400, 500, 600, 700])
        ax.set_yticklabels([])
        # ax.set_title(f"CDS0{WT_num}D_CDS0{HS_num}D\nmannwhitneyu_pvalue_{p_val:.2f}", fontsize=8)
        plt.show()
        # fig.savefig(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_boxPlot_Domain_size_WithoutLable.png', dpi=300, bbox_inches='tight')