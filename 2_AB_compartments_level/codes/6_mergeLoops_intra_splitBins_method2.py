import pathlib as p
import pandas as pd
from intervaltree import Interval, IntervalTree

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns
plt.style.use('ggplot')


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def calculate_accumulated_counts(df):
    # Create a sorted list of all unique midpoints
    midpoints = sorted(set(df['mid1'].tolist() + df['mid2'].tolist()))

    # Create a result DataFrame with intervals and counts initialized to 0
    intervals_df = pd.DataFrame({
        'start': midpoints[:-1],
        'end': midpoints[1:],
        'count': [0] * (len(midpoints) - 1)
    })
    # print(intervals_df)

    # Accumulate counts for each interval
    for _, row in df.iterrows():
        # print(row['count'])
        # Add the count to all intervals that overlap with the current row's interval
        intervals_df['count'] += (
                                         (intervals_df['start'] >= row['mid1']) &
                                         (intervals_df['end'] <= row['mid2'])
                                 ) * row['count']
        # print(intervals_df)
    intervals_df['count'] = intervals_df['count'].astype(int)
    return intervals_df

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '4_annotateLoops'
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    dest_dir = root.parent / 'results' / '6_mergeLoops_intra_splitBins_method2'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 1

    merged_loop_countThresh = 5
    # merged_loop_span_thresh = 500000


    unique_PET_dir = {'01':11754993, '02':12475429, '04':15846412, '05':15597158, '07':15807711, '08':15270275}



    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        print(num)


        # loops_raw = pd.read_csv(src_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone_addType.txt', sep='\t')
        # print(loops_raw.shape)
        # print(loops_raw.head())
        # loops = loops_raw.copy()
        # loops['mid1'] = loops_raw.apply(lambda row: int(0.5*(row['start1']+row['end1'])),axis=1)
        # loops['mid2'] = loops_raw.apply(lambda row: int(0.5 * (row['start2'] + row['end2'])), axis=1)
        # print(loops.shape)
        # print(loops.head())
        # loops.to_csv(dest_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone_addType_addMids.txt', index=None, sep='\t')
        loops = pd.read_csv(dest_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone_addType_addMids.txt', sep='\t')
        # print(loops.shape)
        # print(loops.head(10))

        # print(loops['type'].value_counts())
        loops_filterType = loops[loops['type'] != 'other']
        # # print(loops_filterType.shape)
        # # print(loops_filterType['type'].value_counts())
        loops_filterType_filterMergedCount = loops_filterType[loops_filterType['count'] >= merged_loop_countThresh]
        # print(loops_filterType_filterMergedCount.shape)
        # print(loops_filterType_filterMergedCount.head())
        intra_A = loops_filterType_filterMergedCount[loops_filterType_filterMergedCount['type'] == 'intra_A']
        # print(intra_A.shape)
        # print(intra_A.head())
        intra_A_accumulated = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr4']
        for chro in chro_list:
            print(chro)
            intra_A_filterChro = intra_A[intra_A['chromosome1'] == chro]
            print(intra_A_filterChro.shape)
            print(intra_A_filterChro.head(10))
            intra_A_filterChro_accumulated = calculate_accumulated_counts(intra_A_filterChro)
            intra_A_filterChro_accumulated.insert(0, 'chromosome', chro)
            print(intra_A_filterChro_accumulated.shape)
            print(intra_A_filterChro_accumulated.head(10))
            intra_A_accumulated = pd.concat([intra_A_accumulated, intra_A_filterChro_accumulated])
        print(intra_A_accumulated.shape)
        intra_A_accumulated.to_csv(
            dest_dir / f'CDS0{num}D_intra_A_PETcount_thresh{PETcount_thresh}_merged_loop_countThresh{merged_loop_countThresh}_accumulatedCount.bedgraph',
            index=None, sep='\t')
        
        intra_B = loops_filterType_filterMergedCount[loops_filterType_filterMergedCount['type'] == 'intra_B']
        # print(intra_B.shape)
        # print(intra_B.head())
        intra_B_accumulated = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr4']
        for chro in chro_list:
            print(chro)
            intra_B_filterChro = intra_B[intra_B['chromosome1'] == chro]
            print(intra_B_filterChro.shape)
            print(intra_B_filterChro.head(10))
            intra_B_filterChro_accumulated = calculate_accumulated_counts(intra_B_filterChro)
            intra_B_filterChro_accumulated.insert(0, 'chromosome', chro)
            print(intra_B_filterChro_accumulated.shape)
            print(intra_B_filterChro_accumulated.head(10))
            intra_B_accumulated = pd.concat([intra_B_accumulated, intra_B_filterChro_accumulated])
        print(intra_B_accumulated.shape)
        intra_B_accumulated.to_csv(
            dest_dir / f'CDS0{num}D_intra_B_PETcount_thresh{PETcount_thresh}_merged_loop_countThresh{merged_loop_countThresh}_accumulatedCount.bedgraph',
            index=None, sep='\t')














