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

# def calculate_accumulated_counts(df):
#     # Create a sorted list of all unique midpoints
#     midpoints = sorted(set(df['mid1'].tolist() + df['mid2'].tolist()))
#
#     # Create a result DataFrame with intervals and counts initialized to 0
#     intervals_df = pd.DataFrame({
#         'start': midpoints[:-1],
#         'end': midpoints[1:],
#         'count': [0] * (len(midpoints) - 1),
#         'loop_number': [[] for _ in range(len(midpoints) - 1)]  # Initialize with empty lists
#     })
#
#     # Accumulate counts for each interval
#     for loop_number, row in enumerate(df.itertuples(), 1):  # Using enumerate to get loop number, starting from 1
#         # Find intervals that overlap with the current row's interval
#         overlap_mask = (intervals_df['start'] < row.mid2) & (intervals_df['end'] > row.mid1)
#
#         # Add the count to all intervals that overlap with the current row's interval
#         intervals_df.loc[overlap_mask, 'count'] += row.count
#
#         # Append the loop number to the loop_number column for intervals that overlap
#         for index in intervals_df[overlap_mask].index:
#             intervals_df.at[index, 'loop_number'].append(loop_number)
#
#     return intervals_df


# def calculate_accumulated_counts(df):
#     # Create a sorted list of all unique midpoints
#     midpoints = sorted(set(df['mid1'].tolist() + df['mid2'].tolist()))
#
#     # Create a result DataFrame with intervals and counts initialized to 0
#     intervals_df = pd.DataFrame({
#         'start': midpoints[:-1],
#         'end': midpoints[1:],
#         'count': [0] * (len(midpoints) - 1),
#         'loop_count': [0] * (len(midpoints) - 1)  # Initialize with zeros
#     })
#
#     # Accumulate counts for each interval
#     for loop_number, row in enumerate(df.itertuples(), 1):  # Using enumerate to get loop number, starting from 1
#         # Find intervals that overlap with the current row's interval
#         overlap_mask = (intervals_df['start'] < row.mid2) & (intervals_df['end'] > row.mid1)
#
#         # Add the count to all intervals that overlap with the current row's interval
#         intervals_df.loc[overlap_mask, 'count'] += row.count
#
#         # Increment the loop_count for intervals that overlap
#         intervals_df.loc[overlap_mask, 'loop_count'] += 1
#
#     return intervals_df

# def calculate_accumulated_counts(df):
#     # Create a sorted list of all unique midpoints
#     midpoints = sorted(set(df['mid1'].tolist() + df['mid2'].tolist()))
#
#     # Create a result DataFrame with intervals and counts initialized to 0
#     intervals_df = pd.DataFrame({
#         'start': midpoints[:-1],
#         'end': midpoints[1:],
#         'count': [0] * (len(midpoints) - 1),
#         'loop_count': [0] * (len(midpoints) - 1)  # Initialize with zeros
#     })
#
#     # Accumulate counts for each interval
#     for loop_number, row in enumerate(df.itertuples(), 1):  # Using enumerate to get loop number, starting from 1
#         # Find intervals that overlap with the current row's interval
#         overlap_mask = (intervals_df['start'] < row.mid2) & (intervals_df['end'] > row.mid1)
#
#         # Add the count to all intervals that overlap with the current row's interval
#         intervals_df.loc[overlap_mask, 'count'] += row.count
#
#         # Increment the loop_count for intervals that overlap
#         intervals_df.loc[overlap_mask, 'loop_count'] += 1
#
#     return intervals_df


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir_cis = root.parent / 'data' / 'combined_cis_cluster_fiterChroms_addLoopSpan'
    dest_dir = root.parent / 'results' / '11_call_interaction_domains_STEP1_accumulatedCount'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 5
    loop_span_thresh = 500000


    num_list = ['07', '08','01', '02', '04', '05']
    # num_list = ['07']

    for num in num_list:
        print(num)
        loops_raw = pd.read_csv(src_dir_cis / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan.txt', sep='\t')
        # print(loops_raw.shape)
        # print(loops_raw.head())
        loops_raw_filterPETcount = loops_raw[loops_raw['count'] >= PETcount_thresh]
        # print(loops_raw_filterPETcount.shape)
        # print(loops_raw_filterPETcount.head())
        loops_raw_filterPETcount_filterLoopspan = loops_raw_filterPETcount[loops_raw_filterPETcount['loop_span'] <= loop_span_thresh]
        loops_raw_filterPETcount_filterLoopspan_addMids = loops_raw_filterPETcount_filterLoopspan.copy()
        loops_raw_filterPETcount_filterLoopspan_addMids['mid1'] = loops_raw_filterPETcount_filterLoopspan.apply(lambda row: int(0.5*(row['start1'] + row['end1'])), axis=1)
        loops_raw_filterPETcount_filterLoopspan_addMids['mid2'] = loops_raw_filterPETcount_filterLoopspan.apply(lambda row: int(0.5 * (row['start2'] + row['end2'])), axis=1)
        # print(loops_raw_filterPETcount_filterLoopspan_addMids.shape)
        # print(loops_raw_filterPETcount_filterLoopspan_addMids.head())

        loops = loops_raw_filterPETcount_filterLoopspan_addMids[['chromosome1', 'mid1', 'mid2', 'count']]
        # print(loops.shape)
        # print(loops.head())

        loops_accumulated = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr2L']
        for chro in chro_list:
            print(chro)
            loops_filterChro = loops[loops['chromosome1'] == chro]
            # print(loops_filterChro.shape)
            # print(loops_filterChro.head(10))
            loops_filterChro_accumulated = calculate_accumulated_counts(loops_filterChro)
            loops_filterChro_accumulated.insert(0, 'chromosome', chro)
            print(loops_filterChro_accumulated.shape)
            # print(loops_filterChro_accumulated.head(10))
            loops_accumulated = pd.concat([loops_accumulated, loops_filterChro_accumulated])
        print(loops_accumulated.shape)
        loops_accumulated.to_csv(dest_dir / f'CDS0{num}D_PETcount_thresh{PETcount_thresh}_loop_span_thresh{loop_span_thresh}_accumulatedCount.bedgraph', header=None, index=None, sep='\t')
        loops_accumulated_filterPET = loops_accumulated[loops_accumulated['count']>0]
        print(loops_accumulated_filterPET.shape)
        loops_accumulated_filterPET.to_csv(dest_dir / f'CDS0{num}D_PETcount_thresh{PETcount_thresh}_loop_span_thresh{loop_span_thresh}_accumulatedCount_filterCount0.bedgraph', header=None, index=None, sep='\t')











