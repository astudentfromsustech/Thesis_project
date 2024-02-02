import pathlib as p
import pandas as pd

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

def merge_intervals(df):
    # Sort the DataFrame by the 'start' column
    sorted_df = df.sort_values(by='start')

    # Initialize variables
    merged_intervals = []
    current_start = None
    current_end = None
    current_count = 0

    for _, row in sorted_df.iterrows():
        # If this is the first interval or if the current interval does not overlap
        if current_end is None or row['start'] > current_end:
            # Save the previous interval (if it exists)
            if current_end is not None:
                merged_intervals.append([current_start, current_end, current_count])

            # Start a new interval
            current_start = row['start']
            current_end = row['end']
            current_count = row['accumulated_count']
        else:
            # Merge the current interval with the existing one and sum the counts
            current_end = max(current_end, row['end'])
            current_count += row['accumulated_count']

    # Add the last interval
    merged_intervals.append([current_start, current_end, current_count])

    # Create a new DataFrame
    merged_df = pd.DataFrame(merged_intervals, columns=['start', 'end', 'total_count'])
    return merged_df

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir_cis = root.parent / 'results' / '11_call_interaction_domains_STEP1_accumulatedCount'
    dest_dir = root.parent / 'results' / '11_call_interaction_domains_STEP2_mergeIntervals_filterPercent'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 5
    loop_span_thresh = 500000
    count_percent_thresh = 0.15


    num_list = ['07', '08','01', '02', '04', '05']
    # num_list = ['07']

    for num in num_list:
        print(num)
        intervals = pd.read_csv(src_dir_cis / f'CDS0{num}D_PETcount_thresh{PETcount_thresh}_loop_span_thresh{loop_span_thresh}_accumulatedCount_filterCount0.bedgraph', header=None, sep='\t')
        intervals.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        # print(intervals.shape)
        # print(intervals.head())

        intervals_merged = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr4']
        for chro in chro_list:
            print(chro)
            intervals_Chro = intervals[intervals['chromosome'] == chro]
            print(intervals_Chro.shape)
            print(intervals_Chro.head(10))

            count_thresh = intervals_Chro['accumulated_count'].quantile(count_percent_thresh)
            print(count_thresh)
            intervals_Chro_Count = intervals_Chro[intervals_Chro['accumulated_count'] >= count_thresh]
            # print(intervals_Chro_Count.head(50))
            merged_intervals_Chro_Count = merge_intervals(intervals_Chro_Count)
            merged_intervals_Chro_Count.insert(0, 'chromosome', chro)
            print(merged_intervals_Chro_Count.shape)
            # print(merged_intervals_Chro_Count.head())
            intervals_merged = pd.concat([intervals_merged, merged_intervals_Chro_Count])
        print(intervals_merged.shape)
        # print(intervals_merged)
        intervals_merged.to_csv(dest_dir / f'CDS0{num}D_PETcount_thresh{PETcount_thresh}_loop_span_thresh{loop_span_thresh}_mergedIntervals_filter{count_percent_thresh}P.bed', header=None, index=None, sep='\t')











