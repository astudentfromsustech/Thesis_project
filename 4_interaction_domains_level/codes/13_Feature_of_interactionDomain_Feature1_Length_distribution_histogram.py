import pathlib as p
import pandas as pd

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)



if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '11_call_interaction_domains_STEP3_dropIntervalLen'
    dest_dir = root.parent / 'results' / '13_Feature_of_interactionDomain_Feature1_Length_distribution_histogram'
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
        print(WT_domain_addSize.shape)
        print(WT_domain_addSize.head())

        HS_domain_raw = pd.read_csv(src_dir / f'CDS0{HS_num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
        HS_domain_raw.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        # print(HS_domain_raw.shape)
        # print(HS_domain_raw.head())
        HS_domain = HS_domain_raw[['chromosome', 'start', 'end']]
        # print(HS_domain.shape)
        # print(HS_domain.head())
        HS_domain_addSize = HS_domain.copy()
        HS_domain_addSize['size'] = HS_domain.apply(lambda row: row['end'] - row['start'], axis=1) / 1000
        print(HS_domain_addSize.shape)
        print(HS_domain_addSize.head())

        bin_edges = np.arange(0, 420, 20, dtype=float)
        bin_edges[-1] = np.inf

        hist1, _ = np.histogram(WT_domain_addSize['size'], bins=bin_edges)
        hist2, _ = np.histogram(HS_domain_addSize['size'], bins=bin_edges)
        print('hist1', hist1)
        print('hist2', hist2)
        #
        # Set width for the bars
        width = 5  # adjust width as needed

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
        #
        plt.bar(bar_positions1, hist1, width=width, align='center', color=WT_color)
        plt.bar(bar_positions2, hist2, width=width, align='center', color=HS_color)
        #
        # # plt.xticks(bin_edges[:-1], [f"{edge}-{bin_edges[i + 1]}" for i, edge in enumerate(bin_edges[:-1])], rotation=45)
        ax.set_xticks([0, 100, 200, 300, 400])
        ax.set_xticklabels('')
        # ax.set_xticklabels([20, 120, 220, 320, 420])
        ax.set_yticks([100, 200, 300, 400])
        ax.set_yticklabels([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        # plt.show()
        fig.savefig(dest_dir / f'{WT_num}_{HS_num}_histogram_Domain_size_withoutLabel.png', dpi=300, bbox_inches='tight')