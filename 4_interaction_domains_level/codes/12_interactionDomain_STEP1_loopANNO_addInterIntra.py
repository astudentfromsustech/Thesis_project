import pathlib as p
import pandas as pd

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)


# def annotate_anchors(loops: pd.DataFrame, domains: pd.DataFrame):
#     loops_anno = loops.copy()
#     for idx, row in loops.iterrows():
#         anchor1_intersected_domains = domains[(domains['start'] <= row['end1']) & (domains['end'] >= row['start1'])]
#         print(anchor1_intersected_domains)
#         anchor1_intersected_domains_addOverlappedLen = anchor1_intersected_domains.copy()
#         anchor1_intersected_domains_addOverlappedLen['overlapped_length'] = anchor1_intersected_domains.apply(lambda x: min(row['end1'], x['end']) - max(row['start1'], x['start']), axis=1)
#         print(anchor1_intersected_domains_addOverlappedLen)
#         pass
def annotate_anchors(row, domains: pd.DataFrame):
    anchor1_intersected_domains = domains[(domains['start'] <= row['end1']) & (domains['end'] >= row['start1'])]
    # print(anchor1_intersected_domains)
    if anchor1_intersected_domains.empty:
        anchor1_ID = 'out'
    else:
        anchor1_intersected_domains_addOverlappedLen = anchor1_intersected_domains.copy()
        anchor1_intersected_domains_addOverlappedLen['overlapped_length'] = anchor1_intersected_domains.apply(lambda x: min(row['end1'], x['end']) - max(row['start1'], x['start']), axis=1)
        # print(anchor1_intersected_domains_addOverlappedLen)
        # print(anchor1_intersected_domains_addOverlappedLen['overlapped_length'].idxmax())
        # anchor1_name = anchor1_intersected_domains_addOverlappedLen.loc[anchor1_intersected_domains_addOverlappedLen['overlapped_length'].idxmax()]['name']
        anchor1_ID = anchor1_intersected_domains_addOverlappedLen.loc[anchor1_intersected_domains_addOverlappedLen['overlapped_length'].idxmax()]['ID']

    anchor2_intersected_domains = domains[(domains['start'] <= row['end2']) & (domains['end'] >= row['start2'])]
    # print(anchor2_intersected_domains)
    if anchor2_intersected_domains.empty:
        anchor2_ID = 'out'
    else:
        anchor2_intersected_domains_addOverlappedLen = anchor2_intersected_domains.copy()
        anchor2_intersected_domains_addOverlappedLen['overlapped_length'] = anchor2_intersected_domains.apply(
            lambda x: min(row['end2'], x['end']) - max(row['start2'], x['start']), axis=1)
        # print(anchor2_intersected_domains_addOverlappedLen)
        # print(anchor2_intersected_domains_addOverlappedLen['overlapped_length'].idxmax())
        # anchor2_name = anchor2_intersected_domains_addOverlappedLen.loc[
        #     anchor2_intersected_domains_addOverlappedLen['overlapped_length'].idxmax()]['name']
        anchor2_ID = anchor2_intersected_domains_addOverlappedLen.loc[anchor2_intersected_domains_addOverlappedLen['overlapped_length'].idxmax()]['ID']
    return anchor1_ID, anchor2_ID

def add_type_column(row):
    if row['anchor1_ID'] != 'out' and row['anchor2_ID'] != 'out':
        if row['anchor1_ID'] == row['anchor2_ID']:
            return 'intra'
        else:
            return 'inter'
    else:
        return 'other'


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '12_interactionDomain_STEP1_loopANNO'

    dest_dir = root.parent / 'results' / '12_interactionDomain_STEP1_loopANNO_addInterIntra'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 2
    loop_span_thresh = 500000


    unique_PET_dir = {'01':11754993, '02':12475429, '04':15846412, '05':15597158, '07':15807711, '08':15270275}


    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        print(num)
        loops = pd.read_csv(src_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone.txt', sep='\t')
        # print(loops.shape)
        # print(loops.head(20))
        loops_addType = loops.copy()
        loops_addType['type'] = loops.apply(add_type_column, axis=1)
        print(loops_addType.shape)
        print(loops_addType.head())

        inter = loops_addType[loops_addType['type'] == 'inter']
        # inter = loops_addType[(loops_addType['type'] == 'inter') | (loops_addType['type'] == 'other')]
        print(inter.shape)
        print(inter.head())

        intra = loops_addType[loops_addType['type'] == 'intra']
        print(intra.shape)
        print(intra.head())

        fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            # spine.set_color('none')
            spine.set_linewidth(1)

        plt.scatter(inter['loop_span'], inter['count'], s=1.5, color='#0000FF', alpha=0.7)
        plt.scatter(intra['loop_span'], intra['count'], s=1.5, color='#000033', alpha=0.7)
        plt.xscale('log')
        plt.yscale('log')
        plt.xticks([1000, 10000, 100000, 1000000, 10000000], ['1K', '10K', '100K', '1M', '10M'])
        plt.yticks([2, 10, 100, 1000], ['2', '10', '100', '1K'])
        plt.xlabel('')
        plt.ylabel('')
        plt.savefig(dest_dir / f'CDS0{num}D__filterPETcount{PETcount_thresh}_filterLoopspanNone_inter_intraLoopspan_PETcount_scatterPlot.png', dpi=300)
        # plt.show()
