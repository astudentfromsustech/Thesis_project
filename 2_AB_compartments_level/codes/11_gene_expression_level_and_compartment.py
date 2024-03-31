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

from scipy.stats import pearsonr, spearmanr, gaussian_kde, ttest_ind



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

def add_type_4_column(df):
    df['type_4'] = df['expression_type'].apply(lambda x: 'unregulated' if x in ['up_NS', 'dn_NS'] else x)
    return df

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '10_DEG_and_compartments'
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    dest_dir = root.parent / 'results' / '11_gene_expression_level_and_compartment'
    dest_dir.mkdir(parents=True, exist_ok=True)

    TPM_thresh = 0
    FC_thresh = 1
    FDR_thresh = 0.01
    promoter_extension = 2000
    #
    bin_size = 10000

    #
    DEGs_addID = pd.read_csv(src_dir / f'all_genes_addPromoter{promoter_extension}_addType_TPMthresh{TPM_thresh}_FC{FC_thresh}_FDR{FDR_thresh}_addCompartment.txt', sep='\t')
    # print(DEGs_addID.shape)
    # print(DEGs_addID.head(100))
    DEGs_addID_addAB = add_ab_column(DEGs_addID)
    print(DEGs_addID_addAB.shape)
    print(DEGs_addID_addAB.head())

    DEGs_addID_addAB_addType4 = add_type_4_column(DEGs_addID_addAB)
    print(DEGs_addID_addAB_addType4.shape)
    print(DEGs_addID_addAB_addType4.head())
    statistical_table = pd.crosstab(DEGs_addID_addAB_addType4['type_4'], DEGs_addID_addAB['AB']).reindex(['unexpressed','unregulated', 'up_regulated','dn_regulated'])
    print(statistical_table)

    # all the expressed genes, WT_A, HS_A, HS_A,HS_B  expression level
    expressed = DEGs_addID_addAB_addType4[DEGs_addID_addAB_addType4['type_4'].isin(['unregulated', 'up_regulated', 'dn_regulated'])]
    print(expressed.shape)
    A = expressed[expressed['AB'] == 'A']
    print(A.shape)
    B = expressed[expressed['AB'] == 'B']
    print(B.shape)

    t_stat_pos, p_value_pos = ttest_ind(A['FPKM'], B['FPKM'])

    print(f"T-statistic: {t_stat_pos}")
    print(f"P-value: {p_value_pos}")

    t_stat_neg, p_value_neg = ttest_ind(A['FPKM_RDS026'],  B['FPKM_RDS026'])

    print(f"T-statistic: {t_stat_neg}")
    print(f"P-value: {p_value_neg}")

    df = pd.DataFrame({'WT_A': np.log10(A['FPKM']), 'WT_B':  np.log10(B['FPKM']), 'HS_A':  np.log10(A['FPKM_RDS026']), 'HS_B': np.log10(B['FPKM_RDS026'])})
    # print(df)
    fig, ax = plt.subplots(figsize=(5, 5), dpi=200)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        # spine.set_color('none')
        spine.set_linewidth(1)
    # # #
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    colors_pal = ['#ed8467', '#7b9ef8', '#f1cab6', '#c0d3f5']
    flierprops = dict(marker='o', markerfacecolor='none', markersize=4, markeredgecolor='black', linestyle='none',
                      markeredgewidth=0.6)
    sns.boxplot(data=df, palette=colors_pal, flierprops=flierprops)
    sns.stripplot(data=df, jitter=True, marker='o', alpha=0.7, color='black', s=1)
    # ax.set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
    # ax.set_yticklabels([])
    # ax.set_xticklabels([])

    # ax.set_title(f"{name}_{chro}_WT_HS\nA_pvalue_{p_value_pos:.2f}_B_pvalue_{p_value_neg:.2f}", fontsize=8)
    # plt.savefig(dest_dir / f'{name}_{chro}_WT_HS_A_pvalue_{p_value_pos:.2f}_B_pvalue_{p_value_neg:.2f}_50k_KR_withoutLabel.png', dpi=300, bbox_inches='tight')
    plt.show()


    ##### DEGs and A/B compartment
    # fig, ax = plt.subplots(figsize=(5, 5), dpi=200)
    # ax.set_facecolor('white')
    # for spine in ax.spines.values():
    #     spine.set_color('black')
    #     # spine.set_color('none')
    #     spine.set_linewidth(1)
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    #
    # colors = ['red', 'blue']
    # statistical_table.plot(kind='bar', stacked=True, ax=ax, color=colors, legend=False)
    #
    # # Adding titles and labels
    #
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    # # ax.set_xticks([])
    # # ax.set_xticklabels([])
    # # ax.set_yticklabels([])
    #
    # # Display the plot
    # plt.show()




    ##### draw the pieplot
    # values = [8440, 5745, 1142, 1010]  # unexpressed, unregulated, down-regulated, up-regulated
    #
    # # colors = ['gray', '#B8860B', '#00008B', '#800000']
    # colors = ['gray', '#EE9A4D', 'lightcoral', 'lightblue']
    # fig, ax = plt.subplots(figsize=(5, 5), dpi=200)
    #
    # # Removing the autopct parameter to avoid displaying percentages
    # plt.pie(values, colors=colors, startangle=90)
    # plt.axis('equal')
    #
    # # Display the plot
    # plt.show()

    ##### draw the barplot





















