import pathlib as p
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')

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

    # num_list = [['01', '02', 'green'], ['04', '05', 'red'], ['07', '08','blue']]
    num_list = [['01', '02', 'green']]
    for num_pair in num_list:
        print(num_pair)
        WT_num = num_pair[0]
        HS_num = num_pair[1]
        color = num_pair[2]

        union_intervals =  pd.read_csv(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals_addID_addCounts_addType.txt', sep='\t')
        print(union_intervals.shape)
        print(union_intervals.head())

