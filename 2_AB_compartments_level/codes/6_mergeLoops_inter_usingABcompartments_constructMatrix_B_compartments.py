import pathlib as p
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns
# import typing as t
plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)


def fill_rpm_dataframe(data_df, matrix_df, min_id_number):
    # Copy the data_df to ensure it's not a slice of another DataFrame
    data_df_copy = data_df.copy()

    # Extract and convert anchor IDs to numeric values for the whole DataFrame before looping
    data_df_copy['anchor1_num'] = data_df_copy['anchor1_ID'].str.extract('(\d+)').astype(int)
    data_df_copy['anchor2_num'] = data_df_copy['anchor2_ID'].str.extract('(\d+)').astype(int)

    # Adjust to use numerical indices based on the DataFrame shape or specific logic for index mapping
    for index, row in data_df_copy.iterrows():
        row_index = row['anchor1_num'] - min_id_number
        column_index = row['anchor2_num'] - min_id_number

        # Ensure the indices are within the bounds of matrix_df before assignment
        if 0 <= row_index < len(matrix_df) and 0 <= column_index < len(matrix_df.columns):
            matrix_df.iloc[row_index, column_index] = row['RPM']
            matrix_df.iloc[column_index, row_index] = row['RPM']

    return matrix_df
    

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    src_dir = root.parent / 'results' / '6_mergeLoops_inter_usingABcompartments'
    dest_dir = root.parent / 'results' / '6_mergeLoops_inter_usingABcompartments_constructMatrix'
    dest_dir.mkdir(parents=True, exist_ok=True)


    merged_loop_countThresh_list = [1, 2, 5]
    # merged_loop_span_thresh = 500000
    for merged_loop_countThresh in merged_loop_countThresh_list:
        bin_size = 10000
        chrom = 'chr4'

        intersected_A = pd.read_csv(anno_dir / f'intersected_A_addID_resolution{bin_size}.bed', header=None, sep='\t')
        intersected_A.columns = ['chromosome', 'start', 'end', 'ID']
        # print(intersected_A.shape)
        # print(intersected_A.head())
        intersected_A_chrom = intersected_A[intersected_A['chromosome'] == chrom]
        # print(intersected_A_chrom.shape)
        # print(intersected_A_chrom)

        intersected_B = pd.read_csv(anno_dir / f'intersected_B_addID_resolution{bin_size}.bed', header=None, sep='\t')
        intersected_B.columns = ['chromosome', 'start', 'end', 'ID']
        # print(intersected_B.shape)
        # print(intersected_B.head())
        intersected_B_chrom = intersected_B[intersected_B['chromosome'] == chrom]
        print(intersected_B_chrom.shape)
        print(intersected_B_chrom.head())

        # inter_type = ['inter_B', 'inter_A']

        pair_num = ['01', '02']
        WT_num = pair_num[0]
        HS_num = pair_num[1]
        print(WT_num)
        print(HS_num)
        #
        WT_loops = pd.read_csv(src_dir / f'CDS0{WT_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe', sep='\t')
        WT_inter_B = WT_loops[WT_loops['type'] == 'inter_B']
        # print(WT_inter_B.shape)
        # print(WT_inter_B.head())
        # print(WT_inter_B['type'].value_counts())
        WT_inter_B_chrom = WT_inter_B[(WT_inter_B['chromosome1'] == chrom) & (WT_inter_B['chromosome2'] == chrom)]
        # print(WT_inter_B_chrom.shape)
        # print(WT_inter_B_chrom.head(10))
        # print(WT_inter_B_chrom['type'].value_counts())


        min_id_number = int(intersected_A_chrom['ID'].str.extract('(\d+)').astype(int).min())
        print(min_id_number)
        WT_inter_B_chrom_matrix = pd.DataFrame(0, index=range(intersected_B_chrom.shape[0]), columns=range(intersected_B_chrom.shape[0]))
        # print(WT_inter_B_chrom_matrix.shape)
        # print(WT_inter_B_chrom_matrix)

        WT_inter_B_chrom_matrix_filled = fill_rpm_dataframe(WT_inter_B_chrom, WT_inter_B_chrom_matrix, min_id_number)
        print(WT_inter_B_chrom_matrix_filled.shape)
        print(WT_inter_B_chrom_matrix_filled)

        HS_loops = pd.read_csv(
            src_dir / f'CDS0{HS_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe',
            sep='\t')
        HS_inter_B = HS_loops[HS_loops['type'] == 'inter_B']
        # print(HS_inter_B.shape)
        # print(HS_inter_B.head())
        # print(HS_inter_B['type'].value_counts())
        HS_inter_B_chrom = HS_inter_B[(HS_inter_B['chromosome1'] == chrom) & (HS_inter_B['chromosome2'] == chrom)]
        # print(HS_inter_B_chrom.shape)
        # print(HS_inter_B_chrom.head(10))
        # print(HS_inter_B_chrom['type'].value_counts())

        min_id_number = int(intersected_B_chrom['ID'].str.extract('(\d+)').astype(int).min())
        print(min_id_number)
        HS_inter_B_chrom_matrix = pd.DataFrame(0, index=range(intersected_B_chrom.shape[0]),
                                           columns=range(intersected_B_chrom.shape[0]))
        # print(HS_inter_B_chrom_matrix.shape)
        # print(HS_inter_B_chrom_matrix)

        HS_inter_B_chrom_matrix_filled = fill_rpm_dataframe(HS_inter_B_chrom, HS_inter_B_chrom_matrix, min_id_number)
        print(HS_inter_B_chrom_matrix_filled.shape)
        print(HS_inter_B_chrom_matrix_filled)

        minus = np.log2((HS_inter_B_chrom_matrix_filled + 1) / (WT_inter_B_chrom_matrix_filled + 1))
        print(minus.shape)
        print(minus)

        plt.figure(figsize=(8, 8),dpi=200)  # Set the figure size as desired
        sns.heatmap(minus, annot=False, cmap='bwr', cbar=True, center=0)

        ax = plt.gca()  # Gets the current Axes instance on the current figure
        # ax.xaxis.tick_top()  # Move the X-axis ticks to the top
        # ax.yaxis.tick_right()  # Move the Y-axis ticks to the right
        #
        # # Optionally, if you also want to move the labels to the top and right, you can use:
        # ax.xaxis.set_label_position('top')  # Move the X-axis label to the top
        # ax.yaxis.set_label_position('right')  # Move the Y-axis label to the right
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1)

        # ax.set_xticks(range(0, len(minus.columns), 20))
        # ax.set_yticks(range(0, len(minus.index), 20))
        #
        # # Set the tick labels (optional)
        # ax.set_xticklabels(range(0, len(minus.columns), 20))
        # ax.set_yticklabels(range(0, len(minus.index), 20))

        plt.title(f'log2[(HS_inter_B + 1) / (WT_inter_B + 1)]\n(using all inter_B loops with PET >= {merged_loop_countThresh})')  # Adjust the title position if necessary
        plt.savefig(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_chrom{chrom}_inter_A_loops_with_PETthresh_{merged_loop_countThresh}',dpi=200, bbox_inches='tight')
        plt.show()



















