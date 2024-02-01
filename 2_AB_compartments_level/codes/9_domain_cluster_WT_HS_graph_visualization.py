import pathlib as p
import pandas as pd
from intervaltree import Interval, IntervalTree

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
    edge_dict = {}
    for edge in edges:
        node1, node2, rpm = edge
        if node1 not in edge_dict:
            edge_dict[node1] = {}
        if node2 not in edge_dict:
            edge_dict[node2] = {}
        edge_dict[node1][node2] = (node2, rpm)
        edge_dict[node2][node1] = (node1, rpm)

    def dfs(node, edge_dict, visited, group):
        for neighbor, edge_info in edge_dict[node].items():
            edge = (min(node, neighbor), max(node, neighbor), edge_info[1])
            if edge not in group:
                group.add(edge)
                if neighbor not in visited:
                    visited.add(neighbor)
                    dfs(neighbor, edge_dict, visited, group)

    visited = set()
    groups = []

    for node in edge_dict:
        if node not in visited:
            visited.add(node)
            group = set()
            dfs(node, edge_dict, visited, group)

            # Convert each tuple in the group to a list with string elements
            group_pairs = [[str(item) for item in pair] for pair in group]

            # Sort the pairs and append
            groups.append(sorted(group_pairs, key=lambda x: (x[0], x[1])))

    return groups


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
    dest_dir = root.parent / 'results' / '9_domain_cluster_WT_HS_graph_visualization'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 1
    merged_loop_countThresh = 5
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
    ChIP = stack_and_sort_dataframes(intersected_A, intersected_B)
    # print(ChIP.shape)
    # print(ChIP.head())

    unique_PET_dir = {'01': 11754993, '02': 12475429, '04': 15846412, '05': 15597158, '07': 15807711, '08': 15270275}

    Domain_PETcount_thresh_inter_active = 10
    Domain_distance_thresh_inter_active = 1000000

    # num_list = [['01', '02', '#4D7731', '#98A246'], ['04', '05', '#F50CB8', '#743CA6'],['07', '08', '#1087F4', '#99BDCB']]
    num_list = [['04', '05', '#F50CB8', '#743CA6']]
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
        WT_inter_A_filterDomainCount_filterDomainDistance = WT_inter_A_filterDomainCount[WT_inter_A_filterDomainCount['domain_distance'] <= Domain_distance_thresh_inter_active].head(10)
        print(WT_inter_A_filterDomainCount_filterDomainDistance.shape)
        print(WT_inter_A_filterDomainCount_filterDomainDistance.head(40))

        WT_inter_A_interaction_list = WT_inter_A_filterDomainCount_filterDomainDistance[['anchor1_ID', 'anchor2_ID', 'RPM']].values.tolist()
        print(WT_inter_A_interaction_list)
        G = nx.Graph()
        for edge in WT_inter_A_interaction_list:
            G.add_edge(edge[0], edge[1], weight=edge[2])
        # Define the layout for our nodes
        pos = nx.spring_layout(G)
        # Draw the nodes, edges, and labels
        nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='gray', node_size=2000)
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
        # Show the graph
        plt.title('Graph Plot of Nodes and Interactions')
        plt.show()
        # WT_inter_A_cluster_list = find_connected_groups(WT_inter_A_interaction_list)
        # print(WT_inter_A_cluster_list)
        # WT_inter_A_cluster_df = pd.DataFrame({'cluster': WT_inter_A_cluster_list})
        # print(WT_inter_A_cluster_df)
        # WT_inter_A_cluster_df['num_anchors'], WT_inter_A_cluster_df['num_interaction'], WT_inter_A_cluster_df['cluster_size'] = zip(*WT_inter_A_cluster_df['cluster'].apply(lambda x: calculate_dots_and_edges(x, ChIP)))
        # print(WT_inter_A_cluster_df.shape)
        # print(WT_inter_A_cluster_df.head())
        # WT_inter_A_cluster_df.to_csv(dest_dir / f'CDS0{WT_num}0D_inter_A_cluster_Merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter_active}_Domain_distance{Domain_distance_thresh_inter_active}',
        #                              index=None, sep='\t')

        # HS_inter = pd.read_csv(
        #     src_dir / f'CDS0{HS_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe',
        #     sep='\t')
        # # print(HS_inter.shape)
        # # print(HS_inter.head())
        # HS_inter_A = HS_inter[HS_inter['type'] == 'inter_A']
        # print(HS_inter_A.shape)
        # print(HS_inter_A.head())
        #
        # HS_inter_A_filterDomainCount = HS_inter_A[HS_inter_A['count'] >= Domain_PETcount_thresh_inter_active]
        # # print(HS_inter_A_filterDomainCount.shape)
        # # print(HS_inter_A_filterDomainCount.head())
        # HS_inter_A_filterDomainCount_filterDomainDistance = HS_inter_A_filterDomainCount[
        #     HS_inter_A_filterDomainCount['domain_distance'] <= Domain_distance_thresh_inter_active]
        # # print(HS_inter_A_filterDomainCount_filterDomainDistance.shape)
        # # print(HS_inter_A_filterDomainCount_filterDomainDistance.head(30))
        #
        # HS_inter_A_interaction_list = HS_inter_A_filterDomainCount_filterDomainDistance[
        #     ['anchor1_ID', 'anchor2_ID']].values.tolist()
        # # print(HS_inter_A_interaction_list)
        # HS_inter_A_cluster_list = find_connected_groups(HS_inter_A_interaction_list)
        # # print(HS_inter_A_cluster_list)
        # HS_inter_A_cluster_df = pd.DataFrame({'cluster': HS_inter_A_cluster_list})
        # # print(HS_inter_A_cluster_df)
        # HS_inter_A_cluster_df['num_anchors'], HS_inter_A_cluster_df['num_interaction'], HS_inter_A_cluster_df[
        #     'cluster_size'] = zip(*HS_inter_A_cluster_df['cluster'].apply(lambda x: calculate_dots_and_edges(x, ChIP)))
        # print(HS_inter_A_cluster_df.shape)
        # # print(HS_inter_A_cluster_df)
        # HS_inter_A_cluster_df.to_csv(
        #     dest_dir / f'CDS0{HS_num}0D_inter_A_cluster_Merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter_active}_Domain_distance{Domain_distance_thresh_inter_active}',
        #     index=None, sep='\t')

    # Domain_PETcount_thresh_inter_inactive = 10
    # Domain_distance_thresh_inter_inactive = 500000
    #
    # num_list = [['01', '02', '#4D7731', '#98A246']]
    # for num_pair in num_list:
    #     print(num_pair)
    #     WT_num = num_pair[0]
    #     HS_num = num_pair[1]
    #     WT_color = num_pair[2]
    #     HS_color = num_pair[3]
    #
    #     WT_inter = pd.read_csv(src_dir / f'CDS0{WT_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe', sep='\t')
    #     # print(WT_inter.shape)
    #     # print(WT_inter.head())
    #     WT_inter_B = WT_inter[WT_inter['type'] == 'inter_B']
    #     print(WT_inter_B.shape)
    #     print(WT_inter_B.head())
    #
    #     WT_inter_B_filterDomainCount = WT_inter_B[WT_inter_B['count'] >= Domain_PETcount_thresh_inter_inactive]
    #     # print(WT_inter_B_filterDomainCount.shape)
    #     # print(WT_inter_B_filterDomainCount.head())
    #     WT_inter_B_filterDomainCount_filterDomainDistance = WT_inter_B_filterDomainCount[WT_inter_B_filterDomainCount['domain_distance'] <= Domain_distance_thresh_inter_inactive]
    #     # print(WT_inter_B_filterDomainCount_filterDomainDistance.shape)
    #     # print(WT_inter_B_filterDomainCount_filterDomainDistance.head(30))
    #
    #     WT_inter_B_interaction_list = WT_inter_B_filterDomainCount_filterDomainDistance[['anchor1_ID', 'anchor2_ID']].values.tolist()
    #     # print(WT_inter_B_interaction_list)
    #     WT_inter_B_cluster_list = find_connected_groups(WT_inter_B_interaction_list)
    #     # print(WT_inter_B_cluster_list)
    #     WT_inter_B_cluster_df = pd.DataFrame({'cluster': WT_inter_B_cluster_list})
    #     # print(WT_inter_B_cluster_df)
    #     WT_inter_B_cluster_df['num_anchors'], WT_inter_B_cluster_df['num_interaction'], WT_inter_B_cluster_df['cluster_size'] = zip(*WT_inter_B_cluster_df['cluster'].apply(lambda x: calculate_dots_and_edges(x, ChIP)))
    #     print(WT_inter_B_cluster_df.shape)
    #     # print(WT_inter_B_cluster_df)
    #     WT_inter_B_cluster_df.to_csv(dest_dir / f'CDS0{WT_num}0D_inter_B_cluster_Merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter_inactive}_Domain_distance{Domain_distance_thresh_inter_inactive}',
    #                                  index=None, sep='\t')
    #
    #     HS_inter = pd.read_csv(
    #         src_dir / f'CDS0{HS_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe',
    #         sep='\t')
    #     # print(HS_inter.shape)
    #     # print(HS_inter.head())
    #     HS_inter_B = HS_inter[HS_inter['type'] == 'inter_B']
    #     print(HS_inter_B.shape)
    #     print(HS_inter_B.head())
    #
    #     HS_inter_B_filterDomainCount = HS_inter_B[HS_inter_B['count'] >= Domain_PETcount_thresh_inter_inactive]
    #     # print(HS_inter_B_filterDomainCount.shape)
    #     # print(HS_inter_B_filterDomainCount.head())
    #     HS_inter_B_filterDomainCount_filterDomainDistance = HS_inter_B_filterDomainCount[
    #         HS_inter_B_filterDomainCount['domain_distance'] <= Domain_distance_thresh_inter_inactive]
    #     # print(HS_inter_B_filterDomainCount_filterDomainDistance.shape)
    #     # print(HS_inter_B_filterDomainCount_filterDomainDistance.head(30))
    #
    #     HS_inter_B_interaction_list = HS_inter_B_filterDomainCount_filterDomainDistance[
    #         ['anchor1_ID', 'anchor2_ID']].values.tolist()
    #     # print(HS_inter_B_interaction_list)
    #     HS_inter_B_cluster_list = find_connected_groups(HS_inter_B_interaction_list)
    #     # print(HS_inter_B_cluster_list)
    #     HS_inter_B_cluster_df = pd.DataFrame({'cluster': HS_inter_B_cluster_list})
    #     # print(HS_inter_B_cluster_df)
    #     HS_inter_B_cluster_df['num_anchors'], HS_inter_B_cluster_df['num_interaction'], HS_inter_B_cluster_df[
    #         'cluster_size'] = zip(*HS_inter_B_cluster_df['cluster'].apply(lambda x: calculate_dots_and_edges(x, ChIP)))
    #     print(HS_inter_B_cluster_df.shape)
    #     # print(HS_inter_B_cluster_df)
    #     HS_inter_B_cluster_df.to_csv(
    #         dest_dir / f'CDS0{HS_num}0D_inter_B_cluster_Merged_loop_countThresh{merged_loop_countThresh}_Domain_PETcount{Domain_PETcount_thresh_inter_inactive}_Domain_distance{Domain_distance_thresh_inter_inactive}',
    #         index=None, sep='\t')



















