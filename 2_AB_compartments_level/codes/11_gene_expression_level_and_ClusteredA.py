import pathlib as p
import numpy as np
import pandas as pd
import ast
import seaborn as sns
from itertools import chain

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from matplotlib.colors import LinearSegmentedColormap

from scipy.stats import pearsonr, spearmanr, gaussian_kde, ttest_ind



pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

def get_unique_anchors(cluster):
    seen = set()
    result = []
    for pair in cluster:
        for anchor in pair:
            if anchor not in seen:
                seen.add(anchor)
                result.append(anchor)
    return result

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

def add_type_4_column(df):
    df['type_4'] = df['expression_type'].apply(lambda x: 'unregulated' if x in ['up_NS', 'dn_NS'] else x)
    return df

def get_WT_clustered_ID(compartment_type_ab, WT_clustered_A, WT_clustered_B):
    # Iterate through each row in df2
    if str(compartment_type_ab).startswith('A'):
        for _, row in WT_clustered_A.iterrows():
            # Check if compartment_type_ab is in the unique_anchors list
            if compartment_type_ab in row['unique_anchors']:
                return (row['WT_clustered_A_ID'], row['num_anchors'], row['num_interaction'])  # Return the WT_clustered_A_ID if a match is found
        return (None, None, None)  # Return None if no match is found
    elif str(compartment_type_ab).startswith('B'):
        for _, row in WT_clustered_B.iterrows():
            # Check if compartment_type_ab is in the unique_anchors list
            if compartment_type_ab in row['unique_anchors']:
                return (row['WT_clustered_B_ID'], row['num_anchors'], row['num_interaction'])
    else:
        return (None, None, None)

def get_HS_clustered_ID(compartment_type_ab, HS_clustered_A, HS_clustered_B):
    # Iterate through each row in df2
    if str(compartment_type_ab).startswith('A'):
        for _, row in HS_clustered_A.iterrows():
            # Check if compartment_type_ab is in the unique_anchors list
            if compartment_type_ab in row['unique_anchors']:
                return (row['HS_clustered_A_ID'], row['num_anchors'], row['num_interaction'])  # Return the HS_clustered_A_ID if a match is found
        return (None, None, None)  # Return None if no match is found
    elif str(compartment_type_ab).startswith('B'):
        for _, row in HS_clustered_B.iterrows():
            # Check if compartment_type_ab is in the unique_anchors list
            if compartment_type_ab in row['unique_anchors']:
                return (row['HS_clustered_B_ID'], row['num_anchors'], row['num_interaction'])
    else:
        return (None, None, None)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '10_DEG_and_compartments'
    anno_dir = root.parent / 'results' / '8_domain_cluster_WT_HS'
    dest_dir = root.parent / 'results' / '11_gene_expression_level_and_ClusteredA'
    dest_dir.mkdir(parents=True, exist_ok=True)

    TPM_thresh = 0
    FC_thresh = 1
    FDR_thresh = 0.01
    promoter_extension = 2000
    #
    bin_size = 10000

    Domain_PETcount_thresh = 8

    WT_clustered_A = pd.read_csv(anno_dir / f'CDS0040D_inter_A_cluster_Merged_loop_countThresh5_Domain_PETcount{Domain_PETcount_thresh}_Domain_distance1000000', sep='\t')
    WT_clustered_A['cluster'] = WT_clustered_A['cluster'].apply(ast.literal_eval)
    WT_clustered_A['unique_anchors'] = WT_clustered_A['cluster'].apply(get_unique_anchors)
    WT_clustered_A['WT_clustered_A_ID'] = ['WT_clustered_A_{:03d}'.format(i) for i in range(1, len(WT_clustered_A) + 1)]
    # print(WT_clustered_A.shape)
    # print(WT_clustered_A.head(10))


    HS_clustered_A = pd.read_csv(anno_dir / f'CDS0050D_inter_A_cluster_Merged_loop_countThresh5_Domain_PETcount{Domain_PETcount_thresh}_Domain_distance1000000', sep='\t')
    HS_clustered_A['cluster'] = HS_clustered_A['cluster'].apply(ast.literal_eval)
    HS_clustered_A['unique_anchors'] = HS_clustered_A['cluster'].apply(get_unique_anchors)
    HS_clustered_A['HS_clustered_A_ID'] = ['HS_clustered_A_{:03d}'.format(i) for i in range(1, len(HS_clustered_A) + 1)]
    # print(HS_clustered_A.shape)
    # print(HS_clustered_A.head())

    WT_clustered_B = pd.read_csv(anno_dir / f'CDS0010D_inter_B_cluster_Merged_loop_countThresh5_Domain_PETcount{Domain_PETcount_thresh}_Domain_distance1000000', sep='\t')
    WT_clustered_B['cluster'] = WT_clustered_B['cluster'].apply(ast.literal_eval)
    WT_clustered_B['unique_anchors'] = WT_clustered_B['cluster'].apply(get_unique_anchors)
    WT_clustered_B['WT_clustered_B_ID'] = ['WT_clustered_B_{:03d}'.format(i) for i in range(1, len(WT_clustered_B) + 1)]
    # print(WT_clustered_B.shape)
    # print(WT_clustered_B.head())

    HS_clustered_B = pd.read_csv(anno_dir / f'CDS0020D_inter_B_cluster_Merged_loop_countThresh5_Domain_PETcount{Domain_PETcount_thresh}_Domain_distance1000000', sep='\t')
    HS_clustered_B['cluster'] = HS_clustered_B['cluster'].apply(ast.literal_eval)
    HS_clustered_B['unique_anchors'] = HS_clustered_B['cluster'].apply(get_unique_anchors)
    HS_clustered_B['HS_clustered_B_ID'] = ['HS_clustered_B_{:03d}'.format(i) for i in range(1, len(HS_clustered_B) + 1)]
    # print(HS_clustered_B.shape)
    # print(HS_clustered_B.head())

    DEGs_addID = pd.read_csv(src_dir / f'all_genes_addPromoter{promoter_extension}_addType_TPMthresh{TPM_thresh}_FC{FC_thresh}_FDR{FDR_thresh}_addCompartment.txt', sep='\t')
    # print(DEGs_addID.shape)
    # print(DEGs_addID.head(100))
    DEGs_addID_addAB = add_ab_column(DEGs_addID)
    # print(DEGs_addID_addAB.shape)
    # print(DEGs_addID_addAB.head())

    DEGs_addID_addAB_addType4 = add_type_4_column(DEGs_addID_addAB)
    # print(DEGs_addID_addAB_addType4.shape)
    # print(DEGs_addID_addAB_addType4.head())

    info = DEGs_addID_addAB_addType4[['gene_id', 'name', 'chromosome','FPKM', 'FPKM_RDS026', 'logFC', 'FDR', 'compartment_type', 'AB', 'type_4']]
    # info[['WT_clustered_ID', 'WT_num_anchors', 'WT_num_interaction']] = info['compartment_type'].apply(lambda x: get_WT_clustered_ID(x, WT_clustered_A, WT_clustered_B), result_type='expand')
    expanded_wt = pd.DataFrame(info['compartment_type'].apply(lambda x: get_WT_clustered_ID(x, WT_clustered_A, WT_clustered_B)).tolist(),columns=['WT_clustered_ID', 'WT_num_anchors', 'WT_num_interaction'])
    info = pd.concat([info, expanded_wt], axis=1)
    # info[['HS_clustered_ID', 'HS_num_anchors', 'HS_num_interaction']] = info['compartment_type'].apply(lambda x: get_HS_clustered_ID(x, HS_clustered_A, HS_clustered_B), result_type='expand')
    expanded_hs = pd.DataFrame(info['compartment_type'].apply(lambda x: get_HS_clustered_ID(x, HS_clustered_A, HS_clustered_B)).tolist(), columns=['HS_clustered_ID', 'HS_num_anchors', 'HS_num_interaction'])
    info = pd.concat([info, expanded_hs], axis=1)
    print(info.shape)
    print(info.head(10))

    info.to_csv(dest_dir / f'info_Domain_PETcount{Domain_PETcount_thresh}.txt', index=None, sep='\t')

    # info = pd.read_csv(dest_dir / f'info_Domain_PETcount{Domain_PETcount_thresh}.txt', sep='\t')
    # print(info.shape)
    # print(info.head())
    # print(info['type_4'].value_counts())
    #
    # expressed_list = ['unregulated', 'dn_regulated', 'up_regulated']
    # expressed = info[info['type_4'].isin(expressed_list)]
    # print(expressed.shape)
    # print(expressed.head(3000))


























