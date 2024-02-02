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



if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'combined_cis_cluster_fiterChroms_addLoopSpan'
    interaction_domain_dir = root.parent / 'results' / '11_call_interaction_domains_STEP3_dropIntervalLen'
    dest_dir = root.parent / 'results' / '12_interactionDomain_STEP1_loopANNO'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 2
    # loop_span_thresh = 500000


    unique_PET_dir = {'01':11754993, '02':12475429, '04':15846412, '05':15597158, '07':15807711, '08':15270275}


    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        print(num)
        domains_raw = pd.read_csv(interaction_domain_dir / f'CDS0{num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
        domains_raw.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        domains_raw_addID = domains_raw.copy()
        domains_raw_addID['ID'] = [f'{i:04d}' for i in range(1, len(domains_raw_addID) + 1)]
        # print(domains_raw_addID.shape)
        # print(domains_raw_addID)
        domains = domains_raw_addID[['chromosome', 'start', 'end', 'ID']]
        # print(domains.shape)
        # print(domains.head())

        # loops_raw = pd.read_csv(src_dir / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan.txt', sep='\t')
        # print(loops_raw.shape)
        # print(loops_raw.head())
        # loops_raw_filterPETcount = loops_raw[loops_raw['count'] >= PETcount_thresh]
        # print(loops_raw_filterPETcount.shape)
        # print(loops_raw_filterPETcount.head())
        # # loops_raw_filterPETcount_filterLoopspan = loops_raw_filterPETcount[loops_raw_filterPETcount['loop_span'] <= loop_span_thresh]
        # loops_raw_filterPETcount_filterLoopspan = loops_raw_filterPETcount
        # print(loops_raw_filterPETcount_filterLoopspan.shape)
        # print(loops_raw_filterPETcount_filterLoopspan.head())
        # loops_raw_filterPETcount_filterLoopspan.to_csv(dest_dir / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone.txt',
        #                                                index=None, sep='\t')

        loops = pd.read_csv(dest_dir / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone.txt', sep='\t')
        # print(loops.shape)
        # print(loops.head())

        loops_anno = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr2L']
        for chro in chro_list:
            print(chro)
    #
            loops_chro = loops[loops['chromosome1'] == chro]
            print(loops_chro.shape)
            # print(loops_chro.head())
    #
            domains_chro = domains[domains['chromosome'] == chro]
            print(domains_chro.shape)
            # print(domains_chro.head())
    #
            loops_chro_anno = loops_chro.copy()
            loops_chro_anno[['anchor1_ID', 'anchor2_ID']] = loops_chro.apply(lambda row:annotate_anchors(row, domains_chro), axis=1, result_type='expand')
            print(loops_chro_anno.shape)
            # print(loops_chro_anno)
            loops_anno = pd.concat([loops_anno, loops_chro_anno])
        print(loops_anno.shape)
        # print(loops_anno)
        loops_anno.to_csv(dest_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone.txt', index=None, sep='\t')