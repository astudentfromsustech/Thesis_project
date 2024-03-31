import pathlib as p
import pandas as pd

import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns
import sys
import networkx as nx
# Increase the recursion limit to handle deep recursions
sys.setrecursionlimit(10000)
# import typing as t
plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

def stack_and_sort_dataframes(df1, df2):
    combined_df = pd.concat([df1, df2])
    sorted_df = combined_df.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
    return sorted_df

def add_promoter_regions(df, upstream_length=1000, downstream_length=1000):
    # Initialize new columns with NaN
    df['promoter_start'] = float('nan')
    df['promoter_end'] = float('nan')

    # Calculate promoter regions based on the strand
    for index, row in df.iterrows():
        if row['strand'] == '+':
            df.at[index, 'promoter_start'] = row['start'] - upstream_length
            df.at[index, 'promoter_end'] = row['start'] + downstream_length
        elif row['strand'] == '-':
            df.at[index, 'promoter_start'] = row['end'] - downstream_length
            df.at[index, 'promoter_end'] = row['end'] + upstream_length
    df['promoter_start'] = df['promoter_start'].astype(int)
    df['promoter_end'] = df['promoter_end'].astype(int)
    return df

def add_expressionType(df, TPM_thresh=1, FC_thresh=1, FDR_thresh=0.001):
    df['expression_type'] = None
    for index, row in df.iterrows():
        if pd.isna(row['logFC']) or ((row['TPM'] < TPM_thresh) and (row['TPM_RDS026'] < TPM_thresh)):
            df.at[index, 'expression_type'] = 'unexpressed'
        elif (row['logFC'] >= FC_thresh) and (row['FDR'] <= FDR_thresh):
            df.at[index, 'expression_type'] = 'up_regulated'
        elif ((row['logFC'] >= FC_thresh) and (row['FDR'] > FDR_thresh)) or ((row['logFC'] < FC_thresh) and (row['logFC'] > 0)):
            df.at[index, 'expression_type'] = 'up_NS'
        elif (row['logFC'] <= -FC_thresh) and (row['FDR'] <= FDR_thresh):
            df.at[index, 'expression_type'] = 'dn_regulated'
        elif ((row['logFC'] <= -FC_thresh) and (row['FDR'] > FDR_thresh)) or ((row['logFC'] > -FC_thresh) and (row['logFC'] < 0)):
            df.at[index, 'expression_type'] = 'dn_NS'
    return df

def assign_compartment_type(df1, df2):
    # Create a copy of df1 to avoid modifying the original DataFrame
    df1 = df1.copy()
    # Initialize the new column to None
    df1['compartment_type'] = None
    for index1, row1 in df1.iterrows():
        max_intersect_length = 0
        max_intersect_id = None
        # Filter df2 for matching chromosome
        filtered_df2 = df2[df2['chromosome'] == row1['chromosome']]
        for index2, row2 in filtered_df2.iterrows():
            # Calculate intersection only if chromosomes match
            if row1['chromosome'] == row2['chromosome']:
                intersect_start = max(row1['promoter_start'], row2['start'])
                intersect_end = min(row1['promoter_end'], row2['end'])

                if intersect_start < intersect_end:  # Check if there is an intersection
                    intersect_length = intersect_end - intersect_start
                    if intersect_length > max_intersect_length:
                        max_intersect_length = intersect_length
                        max_intersect_id = row2['ID']
        df1.loc[index1, 'compartment_type'] = max_intersect_id
    return df1

def add_ab_column(df):
    def assign_ab(row):
        if pd.isna(row['compartment_type']):
            return np.nan
        elif row['compartment_type'].startswith('A'):
            return 'A'
        elif row['compartment_type'].startswith('B'):
            return 'B'
        else:
            return np.nan

    df['AB'] = df.apply(assign_ab, axis=1)
    return df

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'annotations'
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    dest_dir = root.parent / 'results' / '10_DEG_and_compartments'
    dest_dir.mkdir(parents=True, exist_ok=True)

    TPM_thresh_list = [0]
    for TPM_thresh in TPM_thresh_list:
        print(TPM_thresh)
        FC_thresh = 1
        FDR_thresh_list = [0.01]
        for FDR_thresh in FDR_thresh_list:
            print(FDR_thresh)
            promoter_extension = 2000
            #
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
            # print(ChIP.shape)
            # print(ChIP.head())
            # ChIP.to_csv(anno_dir / f'stacked_A_B_compartments_addID_resolution10000.bed', index=None, sep='\t')

            all_genes = pd.read_csv(src_dir / f'1_annotated_all_genes_using_dm3_ref_withHighestIsoform_inWT.txt', sep='\t')
            # print(all_genes.shape)
            # print(all_genes.head())
            all_genes_addPromoter = add_promoter_regions(all_genes, upstream_length=promoter_extension, downstream_length=promoter_extension)
            # print(all_genes_addPromoter.shape)
            # print(all_genes_addPromoter.iloc[7900:8100])

            all_genes_addPromoter_addType = add_expressionType(all_genes_addPromoter, TPM_thresh=TPM_thresh, FC_thresh=FC_thresh, FDR_thresh=FDR_thresh)

            all_genes_addPromoter_addType_addCompartment = assign_compartment_type(all_genes_addPromoter_addType, ChIP)
            # print(all_genes_addPromoter_addType_addCompartment.shape)
            # print(all_genes_addPromoter_addType_addCompartment.head(1000))
            all_genes_addPromoter_addType_addCompartment.to_csv(dest_dir / f'all_genes_addPromoter{promoter_extension}_addType_TPMthresh{TPM_thresh}_FC{FC_thresh}_FDR{FDR_thresh}_addCompartment.txt', index=None, sep='\t')

            DEGs_addID = pd.read_csv(dest_dir / f'all_genes_addPromoter{promoter_extension}_addType_TPMthresh{TPM_thresh}_FC{FC_thresh}_FDR{FDR_thresh}_addCompartment.txt', sep='\t')
            # print(DEGs_addID.shape)
            # print(DEGs_addID.head(100))
            DEGs_addID_addAB = add_ab_column(DEGs_addID)
            # print(DEGs_addID_addAB.shape)
            # print(DEGs_addID_addAB.head(100))

            statistical_table = pd.crosstab(DEGs_addID_addAB['expression_type'], DEGs_addID_addAB['AB'])
            print(statistical_table)





















