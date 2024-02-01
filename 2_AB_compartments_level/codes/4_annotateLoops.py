import pathlib as p
import pandas as pd

# import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')



pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def stack_and_sort_dataframes(df1, df2):
    combined_df = pd.concat([df1, df2])
    sorted_df = combined_df.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
    return sorted_df

def annotate_anchors(row, domains: pd.DataFrame):
    # Processing for anchor1
    anchor1_intersected_domains = domains[(domains['start'] <= row['end1']) & (domains['end'] >= row['start1'])]
    if not anchor1_intersected_domains.empty:
        anchor1_intersected_domains['overlapped_length'] = anchor1_intersected_domains.apply(
            lambda x: min(row['end1'], x['end']) - max(row['start1'], x['start']), axis=1)
        anchor1_ID = anchor1_intersected_domains.loc[anchor1_intersected_domains['overlapped_length'].idxmax()]['ID']
    else:
        anchor1_ID = None  # Or some other default value indicating no overlap
    # Processing for anchor2
    anchor2_intersected_domains = domains[(domains['start'] <= row['end2']) & (domains['end'] >= row['start2'])]
    if not anchor2_intersected_domains.empty:
        anchor2_intersected_domains['overlapped_length'] = anchor2_intersected_domains.apply(
            lambda x: min(row['end2'], x['end']) - max(row['start2'], x['start']), axis=1)
        anchor2_ID = anchor2_intersected_domains.loc[anchor2_intersected_domains['overlapped_length'].idxmax()]['ID']
    else:
        anchor2_ID = None  # Or some other default value indicating no overlap
    return anchor1_ID, anchor2_ID

def annotate_type(row):
    anchor1_ID = str(row['anchor1_ID'])
    anchor2_ID = str(row['anchor2_ID'])
    if anchor1_ID.startswith('A') and anchor2_ID.startswith('A'):
        if anchor1_ID == anchor2_ID:
            type = 'intra_A'
        else:
            type = 'inter_A'

    elif anchor1_ID.startswith('B') and anchor2_ID.startswith('B'):
        if anchor1_ID == anchor2_ID:
            type = 'intra_B'
        else:
            type = 'inter_B'

    elif (anchor1_ID.startswith('A') and anchor2_ID.startswith('B')) or (
            anchor1_ID.startswith('B') and anchor2_ID.startswith('A')):
        type = 'inter_mix'

    else:
        type = 'other'
    return type


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '3_combined_cis_cluster_fiterChroms_addLoopSpan'
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    dest_dir = root.parent / 'results' / '4_annotateLoops'
    dest_dir.mkdir(parents=True, exist_ok=True)



    bin_size = 10000


    PETcount_thresh = 1
    # loop_span_thresh = 500000




    intersected_A = pd.read_csv(anno_dir / f'intersected_A_addID_resolution{bin_size}.bed', header=None, sep='\t')
    intersected_A.columns = ['chromosome', 'start', 'end', 'ID']
    # print(intersected_A.shape)
    # print(intersected_A.head())

    intersected_B = pd.read_csv(anno_dir / f'intersected_B_addID_resolution{bin_size}.bed', header=None, sep='\t')
    intersected_B.columns = ['chromosome', 'start', 'end', 'ID']
    # print(intersected_B.shape)
    # print(intersected_B.head())

    ChIP = stack_and_sort_dataframes(intersected_A, intersected_B)
    # print(ChIP.shape)
    # print(ChIP)


    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        print(num)
        loops_raw = pd.read_csv(src_dir / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan.txt', sep='\t')
        print(loops_raw.shape)
        print(loops_raw.head())
        loops_raw_filterPETcount = loops_raw[loops_raw['count'] >= PETcount_thresh]
        print(loops_raw_filterPETcount.shape)
        print(loops_raw_filterPETcount.head())
        # # loops_raw_filterPETcount_filterLoopspan = loops_raw_filterPETcount[loops_raw_filterPETcount['loop_span'] <= loop_span_thresh]
        loops_raw_filterPETcount_filterLoopspan = loops_raw_filterPETcount
        print(loops_raw_filterPETcount_filterLoopspan.shape)
        print(loops_raw_filterPETcount_filterLoopspan.head())
        loops_raw_filterPETcount_filterLoopspan.to_csv(dest_dir / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone.txt',
                                                       index=None, sep='\t')

        loops = pd.read_csv(dest_dir / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone.txt', sep='\t')
        print(loops.shape)
        print(loops.head())
        # #
        loops_anno = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr2L']
        for chro in chro_list:
            print(chro)
        #
            loops_chro = loops[loops['chromosome1'] == chro]
            print(loops_chro.shape)
            print(loops_chro.head())
        #
            ChIP_chro = ChIP[ChIP['chromosome'] == chro]
            print(ChIP_chro.shape)
            print(ChIP_chro.head())
        #
            loops_chro_anno = loops_chro.copy()
            loops_chro_anno[['anchor1_ID', 'anchor2_ID']] = loops_chro.apply(lambda row:annotate_anchors(row, ChIP_chro), axis=1, result_type='expand')
            print(loops_chro_anno.shape)
            print(loops_chro_anno.head())
            loops_anno = pd.concat([loops_anno, loops_chro_anno])
        # print(loops_anno.shape)
        loops_anno.to_csv(dest_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone.txt', index=None, sep='\t')

        loops_anno = pd.read_csv(dest_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone.txt', sep='\t')
        print(loops_anno.shape)
        print(loops_anno.head())
        loops_anno_addType = loops_anno.copy()
        loops_anno_addType['type'] = loops_anno.apply(lambda row: annotate_type(row), axis=1)
        # print(loops_anno_addType.shape)
        # print(loops_anno_addType.tail(1000))
        #
        loops_anno_addType.to_csv(dest_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone_addType.txt', index=None, sep='\t')



















