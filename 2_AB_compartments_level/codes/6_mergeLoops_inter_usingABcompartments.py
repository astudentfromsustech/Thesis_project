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

def stack_and_sort_dataframes(df1, df2):
    combined_df = pd.concat([df1, df2])
    sorted_df = combined_df.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
    return sorted_df

def group_and_sum(df, group_columns, sum_column):
    grouped_df = df.groupby(group_columns)[sum_column].sum().reset_index()
    return grouped_df

def add_domain_info(df1, df2):
    # Merge df1 with df2 for anchor1_ID
    df1 = df1.merge(df2, left_on='anchor1_ID', right_on='ID', how='left')
    df1.rename(columns={'chromosome': 'chromosome1', 'start': 'domain_start1', 'end': 'domain_end1'}, inplace=True)
    df1.drop('ID', axis=1, inplace=True)  # Drop the ID column as it's no longer needed

    # Merge df1 with df2 for anchor2_ID
    df1 = df1.merge(df2, left_on='anchor2_ID', right_on='ID', how='left')
    df1.rename(columns={'chromosome': 'chromosome2', 'start': 'domain_start2', 'end': 'domain_end2'}, inplace=True)
    df1.drop('ID', axis=1, inplace=True)  # Drop the ID column as it's no longer needed
    df1 = df1[['chromosome1', 'domain_start1', 'domain_end1', 'chromosome2', 'domain_start2', 'domain_end2','count', 'anchor1_ID', 'anchor2_ID','type']]
    return df1

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '4_annotateLoops'
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    dest_dir = root.parent / 'results' / '6_mergeLoops_inter_usingABcompartments'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 1

    merged_loop_countThresh = 5
    # merged_loop_span_thresh = 500000


    unique_PET_dir = {'01':11754993, '02':12475429, '04':15846412, '05':15597158, '07':15807711, '08':15270275}

    bin_size = 10000
    intersected_A = pd.read_csv(anno_dir / f'intersected_A_addID_resolution{bin_size}.bed', header=None, sep='\t')
    intersected_A.columns = ['chromosome', 'start', 'end', 'ID']
    # print(intersected_A.shape)
    # print(intersected_A.head())

    intersected_B = pd.read_csv(anno_dir / f'intersected_B_addID_resolution{bin_size}.bed', header=None, sep='\t')
    intersected_B.columns = ['chromosome', 'start', 'end', 'ID']
    # print(intersected_B.shape)
    # print(intersected_B.head())

    ChIP = stack_and_sort_dataframes(intersected_A, intersected_B)
    print(ChIP.shape)
    print(ChIP.head())

    inter_type = ['inter_B', 'inter_A']

    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        print(num)


        loops = pd.read_csv(src_dir / f'ANNO_CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone_addType.txt', sep='\t')
        # print(loops.shape)
        # print(loops.head(10))
        # print(loops['type'].value_counts())
        loops_filterType = loops[loops['type'].isin(inter_type)]
        # print(loops_filterType.shape)
        # print(loops_filterType['type'].value_counts())
        loops_filterType_filterMergedCount = loops_filterType[loops_filterType['count'] >= merged_loop_countThresh]
        print(loops_filterType_filterMergedCount.shape)
        # print(loops_filterType_filterMergedCount.head())
        inter = loops_filterType_filterMergedCount[['count', 'anchor1_ID', 'anchor2_ID', 'type']]
        print(inter.shape)
        # print(inter)

        inter_merged = group_and_sum(inter, ['anchor1_ID', 'anchor2_ID','type'], 'count')
        # print(inter_merged.shape)
        # print(inter_merged)
        inter_merged_addABinfo = add_domain_info(inter_merged, ChIP)
        print(inter_merged_addABinfo.shape)
        print(inter_merged_addABinfo.head())
        # inter_merged_addABinfo.to_csv(dest_dir / f'CDS0{num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}.bedpe', index=None, sep='\t')
        inter_merged_addABinfo_addDistance = inter_merged_addABinfo.copy()
        inter_merged_addABinfo_addDistance['domain_distance'] = inter_merged_addABinfo.apply(lambda row: int(0.5*(row['domain_end2']+row['domain_start2']-row['domain_end1']-row['domain_start1'])), axis=1)
        print(inter_merged_addABinfo_addDistance.shape)
        print(inter_merged_addABinfo_addDistance.head())
        inter_merged_addABinfo_addDistance['RPM'] = round(inter_merged_addABinfo_addDistance['count'] / (unique_PET_dir[num] / 1000000), 2)
        print(inter_merged_addABinfo_addDistance.shape)
        print(inter_merged_addABinfo_addDistance.head())
        inter_merged_addABinfo_addDistance.to_csv(dest_dir / f'CDS0{num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe', index=None, sep='\t')














