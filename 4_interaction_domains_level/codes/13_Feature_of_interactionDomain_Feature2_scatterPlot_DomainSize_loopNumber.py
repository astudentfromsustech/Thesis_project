import pathlib as p
import pandas as pd

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

    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        # loops_raw = pd.read_csv(src_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount2_filterLoopspanNone.txt', sep='\t')
        # print(loops_raw.shape)
        # print(loops_raw.head())
        # loops_raw_addType = loops_raw.copy()
        # loops_raw_addType['type'] = loops_raw.apply(lambda row: annotate_type(row), axis=1)
        # loops_raw_addType.to_csv(src_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount2_filterLoopspanNone_addType.txt', index=None, sep='\t')

        loops_raw_addType = pd.read_csv(src_dir / f'ANNO_CDS001D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount2_filterLoopspanNone_addType.txt', sep='\t')
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

        fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            # spine.set_color('none')
            spine.set_linewidth(1)
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # ax.scatter(np.log10(WT_intra_inactive['loop_span']), np.log10(WT_intra_inactive['count']), s=1.5, color='#000033',alpha=0.7)
        ax.scatter(domains_addSzie['size'] / 1000, domains_addSzie['loop_number'], s=1.5, color='#000033', alpha=0.7)
        # #
        # # # plt.xticks(bin_edges[:-1], [f"{edge}-{bin_edges[i + 1]}" for i, edge in enumerate(bin_edges[:-1])], rotation=45)
        # ax.set_xticks([0, 100, 200, 300, 400])
        # ax.set_xticklabels('')
        # # ax.set_xticklabels([20, 120, 220, 320, 420])
        # ax.set_yticks([100, 200, 300, 400])
        # ax.set_yticklabels([])
        # ax.set_xlabel('')
        # ax.set_ylabel('')
        plt.show()
        # fig.savefig(dest_dir / f'{WT_num}_{HS_num}_histogram_Domain_size_withoutLabel.png', dpi=300, bbox_inches='tight')