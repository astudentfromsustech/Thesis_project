import pathlib as p
import pandas as pd

import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
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

def add_AB_len_column(df1, df2):
    # Merge df1 with df2
    merged_df = pd.merge(df1, df2[['ID', 'AB_len']], left_on='compartment_type', right_on='ID', how='left')
    # Drop the extra 'ID' column from df2
    merged_df.drop('ID', axis=1, inplace=True)
    # Replace missing values in 'AB_len' with NaN
    merged_df['AB_len'].fillna(np.nan, inplace=True)
    return merged_df

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    anno_dir = root.parent / 'results' / '2_index_ABcompartments'
    src_dir = root.parent / 'results' / '10_DEG_and_compartments'
    dest_dir = root.parent / 'results' / '11_gene_expression_level_and_length_of__Acompartments'
    dest_dir.mkdir(parents=True, exist_ok=True)

    bin_size = 10000
    compartments_anno = pd.read_csv(anno_dir / f'stacked_A_B_compartments_addID_resolution{bin_size}.bed', sep='\t')
    # print(compartments_anno.shape)
    # print(compartments_anno.head())
    compartments_anno_addABlen = compartments_anno.copy()
    compartments_anno_addABlen['AB_len'] = compartments_anno.apply(lambda row: row['end'] - row['start'], axis=1)
    # print(compartments_anno_addABlen.shape)
    # print(compartments_anno_addABlen.head())
    # compartments_anno_addABlen.to_csv(anno_dir / f'stacked_A_B_compartments_addID_addABlen_resolution{bin_size}.bed', index=None, sep='\t')

    TPM_thresh = 1
    FC_thresh = 1
    FDR_thresh = 0.001
    promoter_extension = 2000
    DEGs_addID = pd.read_csv(src_dir / f'all_genes_addPromoter{promoter_extension}_addType_TPMthresh{TPM_thresh}_FC{FC_thresh}_FDR{FDR_thresh}_addCompartment.txt', sep='\t')
    DEGs_addID_addAB = add_ab_column(DEGs_addID)
    # print(DEGs_addID_addAB.shape)
    # print(DEGs_addID_addAB.head())
    # statistical_table = pd.crosstab(DEGs_addID_addAB['expression_type'], DEGs_addID_addAB['AB'])
    # print(statistical_table)

    DEGs_addID_addAB_addABlen = add_AB_len_column(DEGs_addID_addAB, compartments_anno_addABlen)
    # print(DEGs_addID_addAB_addABlen.shape)
    # print(DEGs_addID_addAB_addABlen.head(10))

    A_compartments = DEGs_addID_addAB_addABlen[(DEGs_addID_addAB_addABlen['AB'] == 'A')]
    # print(A_compartments.shape)
    # print(A_compartments.head())
    # print(A_compartments['expression_type'].value_counts())

    B_compartments = DEGs_addID_addAB_addABlen[(DEGs_addID_addAB_addABlen['AB'] == 'B')]
    # print(B_compartments.shape)
    # print(B_compartments.head())
    # print(B_compartments['expression_type'].value_counts())

    small_num = 0.001
    A_expressed = A_compartments[A_compartments['expression_type'].isin(['up_NS', 'dn_NS','up_regulated', 'dn_regulated'])]
    A_expressed.loc[A_expressed['FPKM'] == 0, 'FPKM'] += small_num
    A_expressed.loc[A_expressed['FPKM_RDS026'] == 0, 'FPKM_RDS026'] += small_num
    A_expressed.loc[A_expressed['TPM'] == 0, 'TPM'] += small_num
    A_expressed.loc[A_expressed['TPM_RDS026'] == 0, 'TPM_RDS026'] += small_num
    # print(A_expressed.shape)
    # print(A_expressed.head())
    # print(A_expressed['expression_type'].value_counts())

    #
    # A_regulated = A_compartments[A_compartments['expression_type'].isin(['up_regulated', 'dn_regulated'])]
    # print(A_regulated.shape)
    # print(A_regulated.head())
    # print(A_regulated['expression_type'].value_counts())

    A_up_regulated = A_expressed[A_expressed['expression_type'] == 'up_regulated']
    A_dn_regulated = A_expressed[A_expressed['expression_type'] == 'dn_regulated']
    A_up_NS = A_expressed[A_expressed['expression_type'] == 'up_NS']
    A_dn_NS = A_expressed[A_expressed['expression_type'] == 'dn_NS']

    B_expressed = B_compartments[B_compartments['expression_type'].isin(['up_NS', 'dn_NS', 'up_regulated', 'dn_regulated'])]
    B_expressed.loc[B_expressed['FPKM'] == 0, 'FPKM'] += small_num
    B_expressed.loc[B_expressed['FPKM_RDS026'] == 0, 'FPKM_RDS026'] += small_num
    B_expressed.loc[B_expressed['TPM'] == 0, 'TPM'] += small_num
    B_expressed.loc[B_expressed['TPM_RDS026'] == 0, 'TPM_RDS026'] += small_num
    # print(B_expressed.shape)
    # print(B_expressed.head())
    # print(B_expressed['expression_type'].value_counts())
    B_up_regulated = B_expressed[B_expressed['expression_type'] == 'up_regulated']
    B_dn_regulated = B_expressed[B_expressed['expression_type'] == 'dn_regulated']
    B_up_NS = B_expressed[B_expressed['expression_type'] == 'up_NS']
    B_dn_NS = B_expressed[B_expressed['expression_type'] == 'dn_NS']
    # print(B_dn_NS)
    #
    # AB_expressed_data = pd.DataFrame({'A_up_regulated': np.log2(A_up_regulated['FPKM']), 'A_dn_regulated': np.log2(A_dn_regulated['FPKM'])})

    # data = pd.DataFrame({'A_up_regulated': np.log2(A_up_regulated['FPKM']), 'A_dn_regulated': np.log2(A_dn_regulated['FPKM']), 'A_up_NS': np.log2(A_up_NS['FPKM']), 'A_dn_NS': np.log2(A_dn_NS['FPKM']),
    #                      'B_up_regulated': np.log2(B_up_regulated['FPKM']), 'B_dn_regulated': np.log2(B_dn_regulated['FPKM']), 'B_up_NS': np.log2(B_up_NS['FPKM']), 'B_dn_NS': np.log2(B_dn_NS['FPKM'])})

    data = pd.DataFrame({'A_up_regulated': np.log2(A_up_regulated['FPKM']), 'A_up_regulated_HS': np.log2(A_up_regulated['FPKM_RDS026']),
                         'A_dn_regulated': np.log2(A_dn_regulated['FPKM']), 'A_dn_regulated_HS': np.log2(A_dn_regulated['FPKM_RDS026']),
                         'A_up_NS': np.log2(A_up_NS['FPKM']), 'A_up_NS_HS': np.log2(A_up_NS['FPKM_RDS026']),
                         'A_dn_NS': np.log2(A_dn_NS['FPKM']), 'A_dn_NS_HS': np.log2(A_dn_NS['FPKM_RDS026']),
                         'B_up_regulated': np.log2(B_up_regulated['FPKM']), 'B_up_regulated_HS': np.log2(B_up_regulated['FPKM_RDS026']),
                         'B_dn_regulated': np.log2(B_dn_regulated['FPKM']), 'B_dn_regulated_HS': np.log2(B_dn_regulated['FPKM_RDS026']),
                         'B_up_NS': np.log2(B_up_NS['FPKM']), 'B_up_NS_HS': np.log2(B_up_NS['FPKM_RDS026']),
                         'B_dn_NS': np.log2(B_dn_NS['FPKM']), 'B_dn_NS_HS': np.log2(B_dn_NS['FPKM_RDS026'])})

    long_data = data.melt(var_name='Series', value_name='Value')

    fig, ax = plt.subplots(figsize=(10, 6), dpi=200)
    for spine in ax.spines.values():
        spine.set_linewidth(1)
        spine.set_color('black')
    # Set the facecolor of the plot
    ax.set_facecolor('white')
    # sns.violinplot(x='Series', y='Value', data=long_data)
    sns.boxplot(x='Series', y='Value', data=long_data)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=6, rotation=30)
    plt.xlabel('')
    plt.ylabel('FPKM (log2)')
    plt.axhline(y=5, color='r', linestyle='--')
    plt.show()























