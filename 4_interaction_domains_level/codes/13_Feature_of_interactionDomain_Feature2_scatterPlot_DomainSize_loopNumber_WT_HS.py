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

def annotate_type(row):
    if row['anchor1_ID'] != 'out' and row['anchor2_ID'] != 'out':
        if row['anchor1_ID'] == row['anchor2_ID']:
            type = 'intra'
        else:
            type = 'inter'
    else:
        type = 'other'
    return type


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '12_interactionDomain_STEP1_loopANNO'
    domains_dir = root.parent / 'results' / '11_call_interaction_domains_STEP3_dropIntervalLen'
    dest_dir = root.parent / 'results' / '13_Feature_of_interactionDomain_Feature2_scatterPlot_DomainSize_loopNumber'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 5

    # num_list = ['01', '02', '04', '05', '07', '08']
    group_list = [[['01', '#4D7731'], ['02', '#98A246']], [['04', '#F50CB8'], ['05', '#743CA6']], [['07','#1087F4'],  ['08','#99BDCB']]]
    # group_list = [[['04', '#F50CB8'], ['05', '#743CA6']]]
    for group in group_list:

        fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            # spine.set_color('none')
            spine.set_linewidth(1)

        for pair in group:
            print(pair)
            num = pair[0]
            color = pair[1]
            print(num)
            print(color)

            loops_raw_addType = pd.read_csv(src_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount2_filterLoopspanNone_addType.txt', sep='\t')
            # print(loops_raw_addType.shape)
            # print(loops_raw_addType.head())
            loops_raw_addType_filterCount = loops_raw_addType[loops_raw_addType['count'] >= PETcount_thresh]
            # print(loops_raw_addType_filterCount.shape)
            # print(loops_raw_addType_filterCount.head())
            loops_raw_addType_filterCount_filterIntraLoop = loops_raw_addType_filterCount[loops_raw_addType_filterCount['type'] == 'intra']
            # print(loops_raw_addType_filterCount_filterIntraLoop.shape)
            # print(loops_raw_addType_filterCount_filterIntraLoop.head())

            loops = loops_raw_addType_filterCount_filterIntraLoop
            # print(loops)
            grouped = loops.groupby(['chromosome1', 'anchor1_ID']).agg(loop_number=('type', 'size'),  PET_number=('count', 'sum')).reset_index()
            # print(grouped)

            domains = pd.read_csv(domains_dir / f'CDS0{num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
            domains.columns = ['chromosome', 'start', 'end', 'accumulated_count']
            # print(domains.shape)
            # print(domains.head())

            domains['loop_number'] = grouped['loop_number']
            domains['PET_number'] = grouped['PET_number']
            # print(domains.shape)
            # print(domains.head())
            # domains.to_csv(dest_dir / f'CDS0{num}D_filter15P_intervalLen15000_addLoopNumber_addPETnumber.bed', index=None, sep='\t')

            domains_addSzie = domains.copy()
            domains_addSzie['size'] = domains.apply(lambda row: row['end'] - row['start'], axis=1)
            print(domains_addSzie.shape)
            print(domains_addSzie.head())


        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
            ax.scatter(np.log10(domains_addSzie['size'] / 1000), np.log10(domains_addSzie['loop_number']), s=2.5, color=color, alpha=0.7)
        #     ax.scatter((domains_addSzie['size'] / 1000), (domains_addSzie['loop_number']), s=1.5, color=color, alpha=0.7)
        # #
        # fig.savefig(dest_dir / f'{WT_num}_{HS_num}_histogram_Domain_size_withoutLabel.png', dpi=300,
        #             bbox_inches='tight')
        ax.set_xticks([2, 3])
        ax.set_xticklabels([100, 1000])
        ax.set_yticks([1, 2, 3])
        ax.set_yticklabels([10, 100, 1000])
        plt.show()
