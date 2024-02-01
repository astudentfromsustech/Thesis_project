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

def split_into_bins_with_extra(df, num_bins=3, extra_bins=1):
    new_rows = []

    for _, row in df.iterrows():
        start, end = row['start'], row['end']
        interval_size = end - start
        bin_size = interval_size / num_bins

        # Add bins before the interval
        for i in range(extra_bins):
            # print('before', bin_size)
            bin_start = start - (extra_bins - i) * bin_size
            bin_end = bin_start + bin_size
            new_rows.append({
                'chromosome': row['chromosome'],
                'start': int(bin_start),
                'end': int(bin_end),
                'ID': row['ID']
            })

        # Add bins for the interval
        for i in range(num_bins):
            bin_start = start + i * bin_size
            bin_end = bin_start + bin_size
            new_rows.append({
                'chromosome': row['chromosome'],
                'start': int(bin_start),
                'end': int(bin_end),
                'ID': row['ID']
            })

        # Add bins after the interval
        for i in range(extra_bins):
            # print('after', bin_size)
            bin_start = end + i * bin_size
            bin_end = bin_start + bin_size
            new_rows.append({
                'chromosome': row['chromosome'],
                'start': int(bin_start),
                'end': int(bin_end),
                'ID': row['ID']
            })
    return pd.DataFrame(new_rows)

def calculate_accumulated_coverage(df1, df2):
    accumulated_coverages = []

    for _, row1 in df1.iterrows():
        # Initialize the accumulated coverage for the current row of df1
        accumulated_coverage = 0

        # Filter df2 for rows that overlap with the current row of df1
        overlapping_rows = df2[(df2['end'] > row1['start']) &
                               (df2['start'] < row1['end'])]
        # if not overlapping_rows.empty:
        #     print(overlapping_rows.shape)
        #     print(overlapping_rows)
        # Calculate the accumulated coverage
        for _, row2 in overlapping_rows.iterrows():
            # Determine the overlapping interval
            overlap_start = max(row1['start'], row2['start'])
            # print(overlap_start)
            overlap_end = min(row1['end'], row2['end'])
            # print(overlap_end)
            overlap_length = overlap_end - overlap_start
            # print(overlap_length)
            # print(row2['end']-row2['start'])
            # Accumulate the coverage for the overlapping part
            if overlap_length > 0:
                accumulated_coverage += overlap_length * (row2['coverage'] / (row2['end']-row2['start']))

        accumulated_coverages.append(accumulated_coverage)

    # Add the accumulated_coverage column to df1
    df1['accumulated_coverage'] = accumulated_coverages
    df1['bin_size'] = df1.apply(lambda row: row['end'] - row['start'], axis=1)
    return df1
    # return None
def average_elements_by_position(df):
    # Group by the 'ID' column
    grouped = df.groupby('ID')

    # Find the maximum length of groups to know how many averages we need
    max_length = max(grouped.size())

    # Initialize a list of sums and counts, one for each position
    sums = [0] * max_length
    counts = [0] * max_length

    # Iterate through each group
    for _, group in grouped:
        # Iterate over the positions in the group
        for i in range(len(group)):
            sums[i] += group['RPKM'].iloc[i]
            counts[i] += 1

    # Calculate averages
    averages = [sums[i] / counts[i] if counts[i] != 0 else 0 for i in range(max_length)]

    return averages


def moving_average(y_values, window_size):
    # Create an array from the list of y-values
    y = np.array(y_values)

    # Compute the moving average
    cumulative_sum = np.cumsum(y)
    cumulative_sum[window_size:] = cumulative_sum[window_size:] - cumulative_sum[:-window_size]
    smoothed_values = cumulative_sum[window_size - 1:] / window_size

    # Pad the beginning of the smoothed curve with the first value so it has the same length as the original
    smoothed_values = np.concatenate((np.full(window_size - 1, y[0]), smoothed_values))

    return smoothed_values

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '6_mergeLoops_intra_splitBins'
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    dest_dir = root.parent / 'results' / '6_mergeLoops_intra_splitBins_WT_HS_linePlot_method2'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 1

    merged_loop_countThresh = 2
    # merged_loop_span_thresh = 500000

    bin_size = 10000
    intersected_A = pd.read_csv(anno_dir / f'intersected_A_addID_resolution{bin_size}.bed', header=None, sep='\t')
    intersected_A.columns = ['chromosome', 'start', 'end', 'ID']
    # print(intersected_A.shape)
    # print(intersected_A.head())

    intersected_B = pd.read_csv(anno_dir / f'intersected_B_addID_resolution{bin_size}.bed', header=None, sep='\t')
    intersected_B.columns = ['chromosome', 'start', 'end', 'ID']
    # print(intersected_B.shape)
    # print(intersected_B.head())

    num_bins = 300
    extra_bins = 30
    unique_PET_dir = {'01':11754993, '02':12475429, '04':15846412, '05':15597158, '07':15807711, '08':15270275}

    num_list = [['01', '02'], ['04','05'], ['07', '08']]
    # num_list = [['01', '02']]
    for num_pair in num_list:
        num_WT = num_pair[0]
        num_HS = num_pair[1]
        print(num_WT)
        print(num_HS)
    #
        intra_B_accumulated_WT = pd.read_csv(
            src_dir / f'CDS0{num_WT}D_intra_B_PETcount_thresh{PETcount_thresh}_merged_loop_countThresh{merged_loop_countThresh}_accumulatedCount.bedgraph', sep='\t')
        intra_B_accumulated_WT.columns = ['chromosome', 'start', 'end', 'coverage']
        print(intra_B_accumulated_WT.shape)
        print(intra_B_accumulated_WT.head(50))
        #
        intersected_B_allChroms_split_AccumulatedCoverage_addRPKM = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr4']
        for chro in chro_list:
            print(chro)
            intra_B_accumulated_WT_Chro = intra_B_accumulated_WT[intra_B_accumulated_WT['chromosome'] == chro]
            # print(intra_B_accumulated_WT_Chro.shape)
            # print(intra_B_accumulated_WT_Chro.head())

            intersected_B_Chro = intersected_B[intersected_B['chromosome'] == chro]
            # print(intersected_B_Chro.shape)
            # print(intersected_B_Chro)
            intersected_B_Chro_split = split_into_bins_with_extra(intersected_B_Chro, num_bins=num_bins, extra_bins=extra_bins)
            # print(intersected_B_Chro_split.shape)
            # print(intersected_B_Chro_split.head(100))

            intersected_B_Chro_split_AccumulatedCoverage = calculate_accumulated_coverage(intersected_B_Chro_split, intra_B_accumulated_WT_Chro)
            # print(intersected_B_Chro_split_AccumulatedCoverage.shape)
            # print(intersected_B_Chro_split_AccumulatedCoverage.head())
            intersected_B_Chro_split_AccumulatedCoverage_addRPKM = intersected_B_Chro_split_AccumulatedCoverage.copy()
            intersected_B_Chro_split_AccumulatedCoverage_addRPKM['RP'] = intersected_B_Chro_split_AccumulatedCoverage.apply(lambda row: row['accumulated_coverage']/(row['bin_size']/1000), axis=1)
            intersected_B_Chro_split_AccumulatedCoverage_addRPKM['RPKM'] = intersected_B_Chro_split_AccumulatedCoverage.apply(lambda row: row['accumulated_coverage'] / (row['bin_size'] / 1000) / (unique_PET_dir[num_WT]/1000000), axis=1)
            # print(intersected_B_Chro_split_AccumulatedCoverage_addRPKM.shape)
            # print(intersected_B_Chro_split_AccumulatedCoverage_addRPKM.head(60))
            intersected_B_allChroms_split_AccumulatedCoverage_addRPKM = pd.concat([intersected_B_allChroms_split_AccumulatedCoverage_addRPKM, intersected_B_Chro_split_AccumulatedCoverage_addRPKM])
        intersected_B_allChroms_split_AccumulatedCoverage_addRPKM.to_csv(dest_dir / f'CDS0{num_WT}D_intersected_B_allChroms_split_AccumulatedCoverage_addRPKM_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.txt',
                                                                         index=None, sep='\t')
        intersected_B_addRPKM_WT = pd.read_csv(dest_dir / f'CDS0{num_WT}D_intersected_B_allChroms_split_AccumulatedCoverage_addRPKM_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.txt',
                                                                         sep='\t')

        intra_A_accumulated_WT = pd.read_csv(
            src_dir / f'CDS0{num_WT}D_intra_A_PETcount_thresh{PETcount_thresh}_merged_loop_countThresh{merged_loop_countThresh}_accumulatedCount.bedgraph', sep='\t')
        intra_A_accumulated_WT.columns = ['chromosome', 'start', 'end', 'coverage']
        print(intra_A_accumulated_WT.shape)
        print(intra_A_accumulated_WT.head(50))
        #
        intersected_A_allChroms_split_AccumulatedCoverage_addRPKM = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr4']
        for chro in chro_list:
            print(chro)
            intra_A_accumulated_WT_Chro = intra_A_accumulated_WT[intra_A_accumulated_WT['chromosome'] == chro]
            # print(intra_A_accumulated_WT_Chro.shape)
            # print(intra_A_accumulated_WT_Chro.head())

            intersected_A_Chro = intersected_A[intersected_A['chromosome'] == chro]
            # print(intersected_A_Chro.shape)
            # print(intersected_A_Chro)
            intersected_A_Chro_split = split_into_bins_with_extra(intersected_A_Chro, num_bins=num_bins, extra_bins=extra_bins)
            # print(intersected_A_Chro_split.shape)
            # print(intersected_A_Chro_split.head(100))

            intersected_A_Chro_split_AccumulatedCoverage = calculate_accumulated_coverage(intersected_A_Chro_split, intra_A_accumulated_WT_Chro)
            # print(intersected_A_Chro_split_AccumulatedCoverage.shape)
            # print(intersected_A_Chro_split_AccumulatedCoverage.head())
            intersected_A_Chro_split_AccumulatedCoverage_addRPKM = intersected_A_Chro_split_AccumulatedCoverage.copy()
            intersected_A_Chro_split_AccumulatedCoverage_addRPKM['RP'] = intersected_A_Chro_split_AccumulatedCoverage.apply(lambda row: row['accumulated_coverage']/(row['bin_size']/1000), axis=1)
            intersected_A_Chro_split_AccumulatedCoverage_addRPKM['RPKM'] = intersected_A_Chro_split_AccumulatedCoverage.apply(lambda row: row['accumulated_coverage'] / (row['bin_size'] / 1000) / (unique_PET_dir[num_WT]/1000000), axis=1)
            # print(intersected_A_Chro_split_AccumulatedCoverage_addRPKM.shape)
            # print(intersected_A_Chro_split_AccumulatedCoverage_addRPKM.head(60))
            intersected_A_allChroms_split_AccumulatedCoverage_addRPKM = pd.concat([intersected_A_allChroms_split_AccumulatedCoverage_addRPKM, intersected_A_Chro_split_AccumulatedCoverage_addRPKM])
        intersected_A_allChroms_split_AccumulatedCoverage_addRPKM.to_csv(dest_dir / f'CDS0{num_WT}D_intersected_A_allChroms_split_AccumulatedCoverage_addRPKM_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.txt',
                                                                         index=None, sep='\t')
        intersected_A_addRPKM_WT = pd.read_csv(dest_dir / f'CDS0{num_WT}D_intersected_A_allChroms_split_AccumulatedCoverage_addRPKM_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.txt',
                                                                         sep='\t')

        ##########################################################

        intra_B_accumulated_HS = pd.read_csv(
            src_dir / f'CDS0{num_HS}D_intra_B_PETcount_thresh{PETcount_thresh}_merged_loop_countThresh{merged_loop_countThresh}_accumulatedCount.bedgraph',
            sep='\t')
        intra_B_accumulated_HS.columns = ['chromosome', 'start', 'end', 'coverage']
        print(intra_B_accumulated_HS.shape)
        print(intra_B_accumulated_HS.head(50))
        #
        intersected_B_allChroms_split_AccumulatedCoverage_addRPKM = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr4']
        for chro in chro_list:
            print(chro)
            intra_B_accumulated_HS_Chro = intra_B_accumulated_HS[intra_B_accumulated_HS['chromosome'] == chro]
            # print(intra_B_accumulated_HS_Chro.shape)
            # print(intra_B_accumulated_HS_Chro.head())

            intersected_B_Chro = intersected_B[intersected_B['chromosome'] == chro]
            # print(intersected_B_Chro.shape)
            # print(intersected_B_Chro)
            intersected_B_Chro_split = split_into_bins_with_extra(intersected_B_Chro, num_bins=num_bins,
                                                                  extra_bins=extra_bins)
            # print(intersected_B_Chro_split.shape)
            # print(intersected_B_Chro_split.head(100))

            intersected_B_Chro_split_AccumulatedCoverage = calculate_accumulated_coverage(intersected_B_Chro_split,
                                                                                          intra_B_accumulated_HS_Chro)
            # print(intersected_B_Chro_split_AccumulatedCoverage.shape)
            # print(intersected_B_Chro_split_AccumulatedCoverage.head())
            intersected_B_Chro_split_AccumulatedCoverage_addRPKM = intersected_B_Chro_split_AccumulatedCoverage.copy()
            intersected_B_Chro_split_AccumulatedCoverage_addRPKM[
                'RP'] = intersected_B_Chro_split_AccumulatedCoverage.apply(
                lambda row: row['accumulated_coverage'] / (row['bin_size'] / 1000), axis=1)
            intersected_B_Chro_split_AccumulatedCoverage_addRPKM[
                'RPKM'] = intersected_B_Chro_split_AccumulatedCoverage.apply(
                lambda row: row['accumulated_coverage'] / (row['bin_size'] / 1000) / (unique_PET_dir[num_HS] / 1000000),
                axis=1)
            # print(intersected_B_Chro_split_AccumulatedCoverage_addRPKM.shape)
            # print(intersected_B_Chro_split_AccumulatedCoverage_addRPKM.head(60))
            intersected_B_allChroms_split_AccumulatedCoverage_addRPKM = pd.concat(
                [intersected_B_allChroms_split_AccumulatedCoverage_addRPKM,
                 intersected_B_Chro_split_AccumulatedCoverage_addRPKM])
        intersected_B_allChroms_split_AccumulatedCoverage_addRPKM.to_csv(
            dest_dir / f'CDS0{num_HS}D_intersected_B_allChroms_split_AccumulatedCoverage_addRPKM_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.txt',
            index=None, sep='\t')
        intersected_B_addRPKM_HS = pd.read_csv(
            dest_dir / f'CDS0{num_HS}D_intersected_B_allChroms_split_AccumulatedCoverage_addRPKM_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.txt',
            sep='\t')

        intra_A_accumulated_HS = pd.read_csv(
            src_dir / f'CDS0{num_HS}D_intra_A_PETcount_thresh{PETcount_thresh}_merged_loop_countThresh{merged_loop_countThresh}_accumulatedCount.bedgraph',
            sep='\t')
        intra_A_accumulated_HS.columns = ['chromosome', 'start', 'end', 'coverage']
        print(intra_A_accumulated_HS.shape)
        print(intra_A_accumulated_HS.head(50))
        #
        intersected_A_allChroms_split_AccumulatedCoverage_addRPKM = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr4']
        for chro in chro_list:
            print(chro)
            intra_A_accumulated_HS_Chro = intra_A_accumulated_HS[intra_A_accumulated_HS['chromosome'] == chro]
            # print(intra_A_accumulated_HS_Chro.shape)
            # print(intra_A_accumulated_HS_Chro.head())

            intersected_A_Chro = intersected_A[intersected_A['chromosome'] == chro]
            # print(intersected_A_Chro.shape)
            # print(intersected_A_Chro)
            intersected_A_Chro_split = split_into_bins_with_extra(intersected_A_Chro, num_bins=num_bins,
                                                                  extra_bins=extra_bins)
            # print(intersected_A_Chro_split.shape)
            # print(intersected_A_Chro_split.head(100))

            intersected_A_Chro_split_AccumulatedCoverage = calculate_accumulated_coverage(intersected_A_Chro_split,
                                                                                          intra_A_accumulated_HS_Chro)
            # print(intersected_A_Chro_split_AccumulatedCoverage.shape)
            # print(intersected_A_Chro_split_AccumulatedCoverage.head())
            intersected_A_Chro_split_AccumulatedCoverage_addRPKM = intersected_A_Chro_split_AccumulatedCoverage.copy()
            intersected_A_Chro_split_AccumulatedCoverage_addRPKM[
                'RP'] = intersected_A_Chro_split_AccumulatedCoverage.apply(
                lambda row: row['accumulated_coverage'] / (row['bin_size'] / 1000), axis=1)
            intersected_A_Chro_split_AccumulatedCoverage_addRPKM[
                'RPKM'] = intersected_A_Chro_split_AccumulatedCoverage.apply(
                lambda row: row['accumulated_coverage'] / (row['bin_size'] / 1000) / (unique_PET_dir[num_HS] / 1000000),
                axis=1)
            # print(intersected_A_Chro_split_AccumulatedCoverage_addRPKM.shape)
            # print(intersected_A_Chro_split_AccumulatedCoverage_addRPKM.head(60))
            intersected_A_allChroms_split_AccumulatedCoverage_addRPKM = pd.concat(
                [intersected_A_allChroms_split_AccumulatedCoverage_addRPKM,
                 intersected_A_Chro_split_AccumulatedCoverage_addRPKM])
        intersected_A_allChroms_split_AccumulatedCoverage_addRPKM.to_csv(
            dest_dir / f'CDS0{num_HS}D_intersected_A_allChroms_split_AccumulatedCoverage_addRPKM_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.txt',
            index=None, sep='\t')
        intersected_A_addRPKM_HS = pd.read_csv(
            dest_dir / f'CDS0{num_HS}D_intersected_A_allChroms_split_AccumulatedCoverage_addRPKM_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.txt',
            sep='\t')

        # #
        # #
        intersected_B_addRPKM_WT_average = average_elements_by_position(intersected_B_addRPKM_WT)
        # print(pd.Series(intersected_B_addRPKM_WT_average))
        intersected_B_addRPKM_WT_average_smoothed = moving_average(intersected_B_addRPKM_WT_average, 11)

        intersected_B_addRPKM_HS_average = average_elements_by_position(intersected_B_addRPKM_HS)
        intersected_B_addRPKM_HS_average_smoothed = moving_average(intersected_B_addRPKM_HS_average, 11)

        intersected_A_addRPKM_WT_average = average_elements_by_position(intersected_A_addRPKM_WT)
        # print(pd.Series(intersected_A_addRPKM_WT_average))
        intersected_A_addRPKM_WT_average_smoothed = moving_average(intersected_A_addRPKM_WT_average, 11)

        intersected_A_addRPKM_HS_average = average_elements_by_position(intersected_A_addRPKM_HS)
        intersected_A_addRPKM_HS_average_smoothed = moving_average(intersected_A_addRPKM_HS_average, 11)
        # #
        fig, axes = plt.subplots(1,2, figsize=(10, 5), dpi=200)
        for ax in axes:
            ax.set_facecolor('white')
            for spine in ax.spines.values():
                spine.set_color('black')
                spine.set_linewidth(1)
        #
        x_positions = list(range(0, 361, 30))
        print(x_positions)
        x_labels = ['', '0', '30', '60', '90', '120', '150', '180', '210', '240', '270', '300', '']  # Replace with your actual labels
        axes[0].set_xticks(x_positions)
        axes[0].set_xticklabels(x_labels)
        axes[0].plot(intersected_B_addRPKM_WT_average_smoothed, color='black', linewidth=1, label='WT_intra_B')
        axes[0].plot(intersected_B_addRPKM_HS_average_smoothed, color='black', linestyle='--', linewidth=1,label='HS_intra_B')
        axes[0].set_ylabel('RPKM (loop coverage)', fontsize=12)
        axes[0].set_xlabel('Bin', fontsize=12)
        axes[0].legend(loc='upper right', frameon=False)

        axes[1].set_xticks(x_positions)
        axes[1].set_xticklabels(x_labels)
        axes[1].plot(intersected_A_addRPKM_WT_average_smoothed, color='orange', linewidth=1, label='WT_intra_A')
        axes[1].plot(intersected_A_addRPKM_HS_average_smoothed, color='orange', linestyle='--', linewidth=1, label='HS_intra_A')
        axes[1].set_ylabel('RPKM (loop coverage)', fontsize=12)
        axes[1].set_xlabel('Bin', fontsize=12)
        axes[1].legend(loc='upper right', frameon=False)
        # plt.show()
        #
        plt.savefig(dest_dir / f'CDS0{num_WT}D_CDS0{num_HS}D_intersected_B_intersected_A_binNum_{num_bins}_extraNum_{extra_bins}_merged_loop_countThresh{merged_loop_countThresh}.png', dpi=300, bbox_inches='tight')














