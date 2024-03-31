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
    anno_dir = root.parent / 'results' / '8_domain_cluster_WT_HS'
    src_dir = root.parent / 'results' / '10_DEG_and_compartments'
    dest_dir = root.parent / 'results' / '11_DEG_and_length_of_Acompartments_in_a_hub'
    dest_dir.mkdir(parents=True, exist_ok=True)

    TPM_thresh = 1
    FC_thresh = 1
    FDR_thresh = 0.001
    promoter_extension = 2000
    DEGs_addID = pd.read_csv(src_dir / f'all_genes_addPromoter{promoter_extension}_addType_TPMthresh{TPM_thresh}_FC{FC_thresh}_FDR{FDR_thresh}_addCompartment.txt', sep='\t')
    DEGs_addID_addAB = add_ab_column(DEGs_addID)
    print(DEGs_addID_addAB.shape)
    print(DEGs_addID_addAB.head())
    # statistical_table = pd.crosstab(DEGs_addID_addAB['expression_type'], DEGs_addID_addAB['AB'])
    # print(statistical_table)


    merged_loop_countThresh = 5
    Domain_PETcount_thresh_inter = 8
    Domain_distance_thresh_inter = 1000000

    WT_num = '04'
    HS_num = '05'
    WT_A_cluster = pd.read_csv(anno_dir / f'CDS0{WT_num}0D_inter_A_cluster_Merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter}_Domain_distance{Domain_distance_thresh_inter}', sep='\t')
    print(WT_A_cluster.shape)
    print(WT_A_cluster.head())





















