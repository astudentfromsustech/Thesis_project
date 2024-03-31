import pathlib as p
import pandas as pd
from intervaltree import Interval, IntervalTree

import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns

# import typing as t
plt.style.use('ggplot')
# from ggplot import *
# from pandas._libs.tslibs import Timestamp

# import psutil
# from pandarallel import pandarallel
# psutil.cpu_count(logical=False)
# pandarallel.initialize(progress_bar=True, nb_workers=12)

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

def stack_and_sort_dataframes(df1, df2):
    combined_df = pd.concat([df1, df2])
    sorted_df = combined_df.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
    return sorted_df

def find_connected_groups(edges):
    # Helper function to find all connected nodes starting from a specific node
    def dfs(node, group, edge_dict):
        for neighbor in edge_dict[node]:
            if neighbor not in group:
                group.add(neighbor)
                dfs(neighbor, group, edge_dict)

    # Create a dictionary where each node points to its neighbors
    edge_dict = {}
    for edge in edges:
        if edge[0] not in edge_dict:
            edge_dict[edge[0]] = set()
        if edge[1] not in edge_dict:
            edge_dict[edge[1]] = set()
        edge_dict[edge[0]].add(edge[1])
        edge_dict[edge[1]].add(edge[0])

    # Set of all nodes that have already been visited
    visited = set()
    # List of all connected groups
    groups = []

    # Iterate over all nodes
    for node in edge_dict:
        if node not in visited:
            # Start a new group and perform DFS
            group = set()
            dfs(node, group, edge_dict)
            # After DFS, all connected nodes are in 'group'
            groups.append(group)
            # Mark all as visited
            visited.update(group)

    # Convert the sets of groups to lists of lists (pairs)
    group_pairs = []
    for group in groups:
        # Sort the group to maintain consistent order
        group = sorted(group)
        pairs = []
        for node in group:
            # For each node, find all pairs it participates in (node is the smaller end of the pair)
            for connected_node in edge_dict[node]:
                # Since we sorted the group, ensure we only take pairs in the correct order
                if node < connected_node:
                    pairs.append([node, connected_node])
        # Sort pairs and append to the final result
        group_pairs.append(sorted(pairs, key=lambda x: (x[0], x[1])))
    return group_pairs


def calculate_dots_and_edges(cluster, ChIP):

    # Convert each pair (list) in the cluster to a tuple, which is hashable
    edges = [tuple(edge) for edge in cluster]
    # Create a set of dots from the tuples
    dots = set(dot for edge in edges for dot in edge)
    # print(dots)
    # print(max(dots))
    # print(min(dots))
    min_start, min_end = ChIP[ChIP['ID'] == min(dots)].iloc[0]['start'], ChIP[ChIP['ID'] == min(dots)].iloc[0]['end']
    max_start, max_end = ChIP[ChIP['ID'] == max(dots)].iloc[0]['start'], ChIP[ChIP['ID'] == max(dots)].iloc[0]['end']
    # print(min_start, min_end)
    # print(max_start, max_end)
    cluster_size = int(0.5 * (max_end + max_start - min_end - min_start))
    # print(cluster_size)
    # The number of dots is the length of this set
    num_dots = len(dots)
    # The number of edges is the length of the original cluster list
    num_edges = len(cluster)
    return num_dots, num_edges, cluster_size



if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '6_mergeLoops_inter_usingABcompartments'
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    dest_dir = root.parent / 'results' / '8_domain_cluster_WT_HS_forIGV'
    dest_dir.mkdir(parents=True, exist_ok=True)

    # PETcount_thresh = 1
    merged_loop_countThresh = 1
    # merged_loop_span_thresh = 500000


    Domain_PETcount_thresh_inter_active = 0
    Domain_distance_thresh_inter_active = 1000000000

    # num_list = [['01', '02', '#4D7731', '#98A246'], ['04', '05', '#F50CB8', '#743CA6'],['07', '08', '#1087F4', '#99BDCB']]
    num_list = [['04', '05', '#F50CB8', '#743CA6'],['07', '08', '#1087F4', '#99BDCB']]
    for num_pair in num_list:
        print(num_pair)
        WT_num = num_pair[0]
        HS_num = num_pair[1]
        WT_color = num_pair[2]
        HS_color = num_pair[3]

        WT_inter = pd.read_csv(src_dir / f'CDS0{WT_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe', sep='\t')
        # print(WT_inter.shape)
        # print(WT_inter.head())
        WT_inter_A = WT_inter[WT_inter['type'] == 'inter_A']
        print(WT_inter_A.shape)
        print(WT_inter_A.head())

        WT_inter_A_filterDomainCount = WT_inter_A[WT_inter_A['count'] >= Domain_PETcount_thresh_inter_active]
        # print(WT_inter_A_filterDomainCount.shape)
        # print(WT_inter_A_filterDomainCount.head())
        WT_inter_A_filterDomainCount_filterDomainDistance = WT_inter_A_filterDomainCount[WT_inter_A_filterDomainCount['domain_distance'] <= Domain_distance_thresh_inter_active]
        print(WT_inter_A_filterDomainCount_filterDomainDistance.shape)
        print(WT_inter_A_filterDomainCount_filterDomainDistance.head(10))
        WT_inter_A_filterDomainCount_filterDomainDistance.to_csv(
            dest_dir / f'CDS0{WT_num}D_inter_A_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter_active}_Domain_distance{Domain_distance_thresh_inter_active}.bedpe',
            header=None, index=None, sep='\t')


        HS_inter = pd.read_csv(
            src_dir / f'CDS0{HS_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe',
            sep='\t')
        # print(HS_inter.shape)
        # print(HS_inter.head())
        HS_inter_A = HS_inter[HS_inter['type'] == 'inter_A']
        print(HS_inter_A.shape)
        print(HS_inter_A.head())

        HS_inter_A_filterDomainCount = HS_inter_A[HS_inter_A['count'] >= Domain_PETcount_thresh_inter_active]
        # print(HS_inter_A_filterDomainCount.shape)
        # print(HS_inter_A_filterDomainCount.head())
        HS_inter_A_filterDomainCount_filterDomainDistance = HS_inter_A_filterDomainCount[
            HS_inter_A_filterDomainCount['domain_distance'] <= Domain_distance_thresh_inter_active]
        print(HS_inter_A_filterDomainCount_filterDomainDistance.shape)
        print(HS_inter_A_filterDomainCount_filterDomainDistance.head(10))
        HS_inter_A_filterDomainCount_filterDomainDistance.to_csv(
            dest_dir / f'CDS0{HS_num}D_inter_A_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter_active}_Domain_distance{Domain_distance_thresh_inter_active}.bedpe',
            header=None, index=None, sep='\t')

    #
    Domain_PETcount_thresh_inter_inactive = 0
    Domain_distance_thresh_inter_inactive = 1000000000

    num_list = [['01', '02', '#4D7731', '#98A246']]
    for num_pair in num_list:
        print(num_pair)
        WT_num = num_pair[0]
        HS_num = num_pair[1]
        WT_color = num_pair[2]
        HS_color = num_pair[3]

        WT_inter = pd.read_csv(src_dir / f'CDS0{WT_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe', sep='\t')
        # print(WT_inter.shape)
        # print(WT_inter.head())
        WT_inter_B = WT_inter[WT_inter['type'] == 'inter_B']
        print(WT_inter_B.shape)
        print(WT_inter_B.head())

        WT_inter_B_filterDomainCount = WT_inter_B[WT_inter_B['count'] >= Domain_PETcount_thresh_inter_inactive]
        # print(WT_inter_B_filterDomainCount.shape)
        # print(WT_inter_B_filterDomainCount.head())
        WT_inter_B_filterDomainCount_filterDomainDistance = WT_inter_B_filterDomainCount[WT_inter_B_filterDomainCount['domain_distance'] <= Domain_distance_thresh_inter_inactive]
        print(WT_inter_B_filterDomainCount_filterDomainDistance.shape)
        print(WT_inter_B_filterDomainCount_filterDomainDistance.head(10))
        WT_inter_B_filterDomainCount_filterDomainDistance.to_csv(
            dest_dir / f'CDS0{WT_num}D_inter_B_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter_inactive}_Domain_distance{Domain_distance_thresh_inter_inactive}.bedpe',
            header=None, index=None, sep='\t')


        HS_inter = pd.read_csv(
            src_dir / f'CDS0{HS_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe',
            sep='\t')
        # print(HS_inter.shape)
        # print(HS_inter.head())
        HS_inter_B = HS_inter[HS_inter['type'] == 'inter_B']
        print(HS_inter_B.shape)
        print(HS_inter_B.head())

        HS_inter_B_filterDomainCount = HS_inter_B[HS_inter_B['count'] >= Domain_PETcount_thresh_inter_inactive]
        # print(HS_inter_B_filterDomainCount.shape)
        # print(HS_inter_B_filterDomainCount.head())
        HS_inter_B_filterDomainCount_filterDomainDistance = HS_inter_B_filterDomainCount[
            HS_inter_B_filterDomainCount['domain_distance'] <= Domain_distance_thresh_inter_inactive]
        print(HS_inter_B_filterDomainCount_filterDomainDistance.shape)
        print(HS_inter_B_filterDomainCount_filterDomainDistance.head(10))
        HS_inter_B_filterDomainCount_filterDomainDistance.to_csv(
            dest_dir / f'CDS0{HS_num}D_inter_B_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter_inactive}_Domain_distance{Domain_distance_thresh_inter_inactive}.bedpe',
            header=None, index=None, sep='\t')




















