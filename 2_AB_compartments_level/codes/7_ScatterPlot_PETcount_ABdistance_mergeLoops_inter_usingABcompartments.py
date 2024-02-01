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
    src_dir = root.parent / 'results' / '6_mergeLoops_inter_usingABcompartments'
    dest_dir = root.parent / 'results' / '7_ScatterPlot_PETcount_ABdistance_mergeLoops_inter_usingABcompartments'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 1

    merged_loop_countThresh = 5
    # merged_loop_span_thresh = 500000


    unique_PET_dir = {'01':11754993, '02':12475429, '04':15846412, '05':15597158, '07':15807711, '08':15270275}


    num_list = [['01', '02', '#4D7731', '#98A246'], ['04', '05', '#F50CB8', '#743CA6'], ['07', '08', '#1087F4', '#99BDCB']]
    # num_list = [['01', '02','#4D7731','#98A246']]
    for num_pair in num_list:
        print(num_pair)
        WT_num = num_pair[0]
        HS_num = num_pair[1]
        WT_color = num_pair[2]
        HS_color = num_pair[3]

        WT_inter = pd.read_csv(src_dir / f'CDS0{WT_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe', sep='\t')
        print(WT_inter.shape)
        print(WT_inter.head())

        WT_inter_A = WT_inter[WT_inter['type'] == 'inter_A']
        print(WT_inter_A.shape)
        print(WT_inter_A.head())

        WT_inter_B = WT_inter[WT_inter['type'] == 'inter_B']
        print(WT_inter_B.shape)
        print(WT_inter_B.head())

        HS_inter = pd.read_csv(src_dir / f'CDS0{HS_num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe', sep='\t')
        print(HS_inter.shape)
        print(HS_inter.head())

        HS_inter_A = HS_inter[HS_inter['type'] == 'inter_A']
        print(HS_inter_A.shape)
        print(HS_inter_A.head())

        HS_inter_B = HS_inter[HS_inter['type'] == 'inter_B']
        print(HS_inter_B.shape)
        print(HS_inter_B.head())

        plt.rcParams['figure.dpi'] = 300

        # Initialize the JointGrid
        g = sns.JointGrid(x=np.log10(WT_inter_A['domain_distance']),
                          y=np.log10(WT_inter_A['count']),
                          height=5)

        # Plot the WT scatter plot on the JointGrid
        g.plot_joint(plt.scatter,
                     color=WT_color,
                     s=5,
                     edgecolor='white',
                     label='WT inter_A',
                     alpha=0.7)

        # Overlay the HS scatter plot
        g.ax_joint.scatter(np.log10(HS_inter_A['domain_distance']),
                           np.log10(HS_inter_A['count']),
                           color=HS_color,
                           s=5,
                           edgecolor='white',
                           alpha=0.7,
                           label='HS inter_A')

        # Plot the WT density plots on the margins
        sns.kdeplot(x=np.log10(WT_inter_A['domain_distance']),
                    ax=g.ax_marg_x,
                    color=WT_color,
                    linewidth=1)
        sns.kdeplot(y=np.log10(WT_inter_A['count']),
                    ax=g.ax_marg_y,
                    color=WT_color,
                    linewidth=1)

        # Overlay the HS density plots on the margins
        sns.kdeplot(x=np.log10(HS_inter_A['domain_distance']),
                    ax=g.ax_marg_x,
                    color=HS_color,
                    linewidth=1)
        sns.kdeplot(y=np.log10(HS_inter_A['count']),
                    ax=g.ax_marg_y,
                    color=HS_color,
                    linewidth=1)

        # Set the background color of the whole figure and the individual plots
        g.fig.set_facecolor('white')
        g.ax_joint.set_facecolor('white')
        g.ax_joint.set_yticks([1, 2, 3])
        g.ax_joint.set_yticklabels([10, 100, 1000])
        g.ax_joint.set_ylabel('PET count')
        g.ax_joint.set_xticks([4, 5, 6, 7])
        g.ax_joint.set_xticks([4, 5, 6, 7])
        g.ax_joint.set_xticklabels(['10K', '100K', '1M', '10M'])
        g.ax_marg_x.set_facecolor('white')
        g.ax_marg_y.set_facecolor('white')

        # Set the spines to be visible
        for ax in [g.ax_joint, g.ax_marg_x, g.ax_marg_y]:
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color('black')
                spine.set_linewidth(1)

        # Adjust the figure to prevent cut-off titles or labels
        g.ax_joint.legend()
        g.ax_joint.legend(frameon=False, facecolor='none')
        plt.tight_layout()

        # Display the plot
        # plt.show()
        g.fig.savefig(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_inter_A_domainDistance_densityPlot_merged_loop_countThresh{merged_loop_countThresh}.png', dpi=300, bbox_inches='tight')













