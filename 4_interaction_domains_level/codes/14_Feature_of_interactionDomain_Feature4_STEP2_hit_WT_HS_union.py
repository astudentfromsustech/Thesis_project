import pathlib as p
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from matplotlib.colors import LinearSegmentedColormap

from scipy.stats import pearsonr, spearmanr, gaussian_kde, ttest_ind, mannwhitneyu, wilcoxon, ks_2samp

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def union_intervals(df1, df2):
    # Concatenate the two dataframes
    df = pd.concat([df1[['start', 'end']], df2[['start', 'end']]])

    # Sort by the 'start' column
    df = df.sort_values(by='start')

    # Iterate and merge intervals
    merged = []
    for _, row in df.iterrows():
        if not merged or merged[-1][1] < row['start']:
            # If no overlap, add the interval as is
            merged.append([row['start'], row['end']])
        else:
            # If there is an overlap, merge with the last interval
            merged[-1][1] = max(merged[-1][1], row['end'])

    # Convert to DataFrame
    return pd.DataFrame(merged, columns=['start', 'end'])

def assign_id(df1, df2):
    # Create a new column for IDs in df1, initialize with None
    df1['ID'] = None

    # Iterate over each row in df1
    for i, row1 in df1.iterrows():
        # Check each row in df2
        for _, row2 in df2.iterrows():
            # Check if chromosomes match and df1 interval is within df2 interval
            if row1['chromosome'] == row2['chromosome'] and row1['start'] >= row2['start'] and row1['end'] <= row2['end']:
                # Assign the ID from df2 to df1
                df1.at[i, 'ID'] = row2['ID']
                break  # Stop checking once the first match is found
    return df1

def anno_type_WT(WT, df):
    WT['type'] = None
    for i, row1 in WT.iterrows():
        # print(row1['ID'])
        # print(df['ID'].isin([row1['ID']]).any())
        if df['ID'].isin([row1['ID']]).any():
            WT.at[i, 'type'] = 'shared'
        else:
            WT.at[i, 'type'] = 'WT_only'
    return WT

def anno_type_HS(HS, df):
    HS['type'] = None
    for i, row1 in HS.iterrows():
        # print(row1['ID'])
        # print(df['ID'].isin([row1['ID']]).any())
        if df['ID'].isin([row1['ID']]).any():
            HS.at[i, 'type'] = 'shared'
        else:
            HS.at[i, 'type'] = 'HS_only'
    return HS

def count_WT(df1, df2, df3):
    # Count the occurrences of each ID in df2
    id_counts_WT = df2['ID'].value_counts()
    id_counts_HS = df3['ID'].value_counts()
    # Map these counts to df1 based on the ID column
    
    df1['WT_count'] = df1['ID'].map(id_counts_WT)
    df1['WT_count'].fillna(0, inplace=True)
    df1['HS_count'] = df1['ID'].map(id_counts_HS)
    df1['HS_count'].fillna(0, inplace=True)

    return df1
def annotate_merge_in_union(union):
    type_list = []
    for idx, row in union.iterrows():
        if row['WT_count'] == 1 and row['HS_count'] == 1:
            type = 'keep'
        elif row['WT_count'] == 1 and row['HS_count'] == 0:
            type = 'disappear'
        elif row['WT_count'] == 0 and row['HS_count'] == 1:
            type = 'appear'
        elif row['WT_count'] > row['HS_count']:
            type = 'merge'
        else:
            type = 'split'
        type_list.append(type)
    union['type'] = type_list
    return union

def chip_anno(df1):
    ChIP = pd.read_csv(anno_dir / f'active.minji.txt', sep='\t')
    ChIP_addSize = ChIP.copy()
    ChIP_addSize['size'] = ChIP_addSize.apply(lambda row: row['End'] - row['Start'], axis=1)
    # print(ChIP_addSize.shape)
    # print(ChIP_addSize.head())
    active = ChIP_addSize[ChIP_addSize['name'] == 'MJAD3']
    # print(sum(active['size']))
    inactive = ChIP_addSize[ChIP_addSize['name'] == 'MJID4']
    # print(sum(inactive['size']))
    # print(sum(ChIP_addSize['size']))
    ratio = sum(inactive['size']) / sum(active['size'])

    df1 = df1.copy()

    # Initialize the new column
    df1['chip'] = 'inactive'  # Default value, can be changed later

    # Iterate over each row in df1
    for index, row in df1.iterrows():
        chromosome1, start1, end1 = row['chromosome'], row['start'], row['end']
        mjad3_length, mjid4_length = 0, 0

        # Filter ChIP_addSize for rows with the same chromosome
        ChIP_addSize_filtered = ChIP_addSize[ChIP_addSize['#Chromosome'] == chromosome1]

        # Check against all rows in filtered ChIP_addSize
        for _, row2 in ChIP_addSize_filtered.iterrows():
            start2, end2 = row2['Start'], row2['End']
            # Check for intersection
            if start1 < end2 and start2 < end1:
                intersect_len = min(end1, end2) - max(start1, start2)
                # Accumulate the lengths
                if row2['name'] == 'MJAD3':
                    mjad3_length += intersect_len
                elif row2['name'] == 'MJID4':
                    mjid4_length += intersect_len

        # Assign 'active' or 'inactive'
        # print(mjad3_length)
        # print(mjid4_length / ratio)
        df1.loc[index, 'chip'] = 'active' if mjad3_length > (mjid4_length / ratio) else 'inactive'

    return df1

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    anno_dir = root.parent / 'annotations'
    src_dir = root.parent / 'results' / '11_call_interaction_domains_STEP3_dropIntervalLen'
    dest_dir = root.parent / 'results' / '14_Feature_of_interactionDomain_Feature4_STEP1_WT_HS_union'
    dest_dir.mkdir(parents=True, exist_ok=True)

    num_list = [['01', '02', 'green'], ['04', '05', 'red'], ['07', '08','blue']]
    # num_list = [['01', '02', 'green']]
    for num_pair in num_list:
        print(num_pair)
        WT_num = num_pair[0]
        HS_num = num_pair[1]
        color = num_pair[2]


        WT_domain_raw = pd.read_csv(src_dir / f'CDS0{WT_num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
        WT_domain_raw.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        # print(WT_domain_raw.shape)
        # print(WT_domain_raw.head())
        WT_domain = WT_domain_raw[['chromosome', 'start', 'end']]
        # print(WT_domain.shape)
        # print(WT_domain.head())


        HS_domain_raw = pd.read_csv(src_dir / f'CDS0{HS_num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
        HS_domain_raw.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        # print(HS_domain_raw.shape)
        # print(HS_domain_raw.head())
        HS_domain = HS_domain_raw[['chromosome', 'start', 'end']]
        # print(HS_domain.shape)
        # print(HS_domain.head())

        # union = pd.DataFrame()
        # chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # # chro_list = ['chr2L']
        # for chro in chro_list:
        #     print(chro)
        #     WT_domain_chro = WT_domain[WT_domain['chromosome'] == chro]
        #     # print(WT_domain_chro.shape)
        #     # print(WT_domain_chro.head())
        #
        #     HS_domain_chro = HS_domain[HS_domain['chromosome'] == chro]
        #     # print(HS_domain_chro.shape)
        #     # print(HS_domain_chro.head())
        #
        #     union_chro = union_intervals(WT_domain_chro, HS_domain_chro)
        #     print(union_chro.shape)
        #     union_chro.insert(0, 'chromosome', chro)
        #     union = pd.concat([union, union_chro])
        # print(union.shape)
        # union['ID'] = [f'{i:03d}' for i in range(1, len(union) + 1)]
        # print(union.shape)
        # print(union)
        # union.to_csv(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals.txt', index=None, sep='\t')
        union = pd.read_csv(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals.txt', dtype={'ID': str}, sep='\t')
        # print(union.shape)
        # print(union.head())

        # WT_domain_addAnno = assign_id(WT_domain, union)
        # WT_domain_addAnno.to_csv(dest_dir / f'CDS0{WT_num}D_hit_union_anno.txt', index=None, sep='\t')
        WT_domain_addAnno = pd.read_csv(dest_dir / f'CDS0{WT_num}D_hit_union_anno.txt', dtype={'ID': str}, sep='\t')
        # print(WT_domain_addAnno.shape)
        # print(WT_domain_addAnno.head())

        # HS_domain_addAnno = assign_id(HS_domain, union)
        # HS_domain_addAnno.to_csv(dest_dir / f'CDS0{HS_num}D_hit_union_anno.txt', index=None, sep='\t')
        HS_domain_addAnno = pd.read_csv(dest_dir / f'CDS0{HS_num}D_hit_union_anno.txt', dtype={'ID': str}, sep='\t')
        # print(HS_domain_addAnno.shape)
        # print(HS_domain_addAnno.head(20))


                        # # # WT_domain_addAnno_addType = anno_type_WT(WT_domain_addAnno, HS_domain_addAnno)
                        # # # WT_domain_addAnno_addType.to_csv(dest_dir / f'CDS0{WT_num}D_hit_union_anno_addType.txt', index=None, sep='\t')
                        # # # # print(WT_domain_addAnno_addType)
                        # # # HS_domain_addAnno_addType = anno_type_HS(HS_domain_addAnno, WT_domain_addAnno)
                        # # # HS_domain_addAnno_addType.to_csv(dest_dir / f'CDS0{HS_num}D_hit_union_anno_addType.txt', index=None, sep='\t')
                        # # # # print(HS_domain_addAnno_addType)

        union_addCount = count_WT(union, WT_domain_addAnno, HS_domain_addAnno)
        # print(union_addCount)

        union_addCount_addType = annotate_merge_in_union(union_addCount)
        # print(union_addCount_addType.shape)
        # print(union_addCount_addType.head())
        # union_addCount_addType.to_csv(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals_addID_addCounts_addType.txt', index=None, sep='\t')

        #
                # # # fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
                # # # ax.set_facecolor('white')
                # # # for spine in ax.spines.values():
                # # #     spine.set_color('black')
                # # #     # spine.set_color('none')
                # # #     spine.set_linewidth(1)
                # # # ax.spines['top'].set_visible(False)
                # # # ax.spines['right'].set_visible(False)
                # # #
                # # # type_counts = union_addCount_addType['type'].value_counts().reindex(['keep', 'appear', 'disappear', 'merge', 'split'])
                # # # print(sum(type_counts))
                # # # print(type_counts)
                # # # # type_counts.plot(kind='bar', color=color)
                # # # # # plt.title('Counts of Types')
                # # # # plt.xlabel('')
                # # # # plt.ylabel('Count')
                # # # # plt.xticks(rotation=25)
                # # # # # plt.show()
                # # # # plt.savefig(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals_addID_addCounts_addType_barPlot.png', dpi=300, bbox_inches='tight')


        # union_addCount_addType_addChip = chip_anno(union_addCount_addType)
        # union_addCount_addType_addChip.to_csv(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals_addID_addCounts_addType_addChip.txt', index=None, sep='\t')

        union_addCount_addType_addChip = pd.read_csv(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals_addID_addCounts_addType_addChip.txt', sep='\t')
        print(union_addCount_addType_addChip.shape)
        # print(union_addCount_addType_addChip)
        stat_table = pd.crosstab(union_addCount_addType_addChip['type'], union_addCount_addType_addChip['chip']).reindex(['keep', 'appear', 'disappear', 'merge', 'split'])
        print(stat_table)

        colors = ['orange', 'black']
        fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            # spine.set_color('none')
            spine.set_linewidth(1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)


        stat_table.plot(kind='bar', ax=ax, color=colors)
        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('Count')
        ax.set_xticklabels(stat_table.index, rotation=25)
        ax.get_legend().remove()
        # plt.show()
        plt.savefig(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals_addID_addCounts_addType_addChip_barPlot.png',dpi=300, bbox_inches='tight')