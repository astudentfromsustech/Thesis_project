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

from scipy.stats import pearsonr, spearmanr, gaussian_kde, wilcoxon, ttest_rel

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

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

# from scipy.stats import pearsonr, spearmanr, gaussian_kde

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def get_expression_data(row, df2):
    if row['expression'] == 'unexpressed':
        return None, None, None, None, None, None
    else:
        match = df2[df2['gene_id'] == row['gene_symbol']]
        if not match.empty:
            return match.iloc[0]['TPM'], match.iloc[0]['TPM_RDS026'], match.iloc[0]['FPKM'], match.iloc[0]['FPKM_RDS026'], match.iloc[0]['logFC'], match.iloc[0]['FDR']
        else:
            return None, None, None, None, None, None

def get_domain_ID_type(row, df2):
    same_chromosome = df2[df2['chromosome'] == row['chromosome']]
    # Check for intersection
    for _, row2 in same_chromosome.iterrows():
        if row['P_left'] <= row2['end'] and row['P_right'] >= row2['start']:
            return row2['ID'], row2['start'], row2['end'], row2['type']
    return None,None,None,None

def map_symbols_to_ids(series, dataframe):
    # Create a mapping dictionary from the dataframe
    mapping_dict = dataframe.set_index('symbol')['ID'].to_dict()

    # Map the series values to the corresponding IDs
    mapped_series = series.map(mapping_dict)

    # Drop the elements where no matching ID is found
    mapped_series.dropna(inplace=True)

    return mapped_series

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data'
    domain_dir = root.parent / 'results' / '14_Feature_of_interactionDomain_Feature4_STEP1_WT_HS_union'
    anno_dir = root.parent / 'annotations'
    dest_dir = root.parent / 'results' / '22_RNA-seq_with_largestExpression_isoform_add5Category_GO'
    dest_dir.mkdir(exist_ok=True, parents=True)

    TPM_thresh = 0
    logFC_thresh = 1
    FDR_thresh = 0.001
    IDs = pd.read_csv(anno_dir / f'final_dm3_Ensembl_ID_and_geneSymbol.txt', header=None, sep='\t')
    IDs.columns = ['ID', 'symbol']
    print(IDs.shape)
    print(IDs.head())

    # UCSC_ref_addLabel_addType_Chrolist = pd.read_csv(src_dir / 'RNA-seq' / f'UCSC_ref_addLabel_addType_chroList_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}.txt', sep='\t')
    # print(UCSC_ref_addLabel_addType_Chrolist.shape)
    # print(UCSC_ref_addLabel_addType_Chrolist.head())
    # # # #
    # # # # expressed = UCSC_ref_addLabel_addType[UCSC_ref_addLabel_addType['expression'] == 'expressed']
    # # # # print(expressed.shape)
    # # # # print(expressed.head())
    # # # # up_DEGs = UCSC_ref_addLabel_addType[UCSC_ref_addLabel_addType['type'] == 'up_regulated']
    # # # # print(up_DEGs.shape)
    # # # # print(up_DEGs.head())
    # # # # down_DEGs = UCSC_ref_addLabel_addType[UCSC_ref_addLabel_addType['type'] == 'down_regulated']
    # # # # print(down_DEGs.shape)
    # # # # print(down_DEGs.head())

    num_list = [['01', '02', '#4D7731', '#98A246'], ['04', '05', '#F50CB8', '#743CA6'], ['07', '08', '#1087F4', '#99BDCB']]
    # num_list = [['01', '02', '#4D7731', '#98A246']]
    for num_pair in num_list:
        print(num_pair)
        WT_num = num_pair[0]
        HS_num = num_pair[1]
        WT_color = num_pair[2]
        HS_color = num_pair[3]

    #     union_intervals =  pd.read_csv(domain_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals_addID_addCounts_addType.txt', sep='\t', dtype={'ID': str})
    #     print(union_intervals.shape)
    #     print(union_intervals.head())
    # #
    #     UCSC_ref_addLabel_addType_Chrolist_addCategory = UCSC_ref_addLabel_addType_Chrolist.copy()
    #     UCSC_ref_addLabel_addType_Chrolist_addCategory[['ID', 'domain_start','domain_end','category']] = UCSC_ref_addLabel_addType_Chrolist.apply(lambda row: get_domain_ID_type(row, union_intervals), axis=1, result_type='expand')
    #     print(UCSC_ref_addLabel_addType_Chrolist_addCategory.shape)
    #     print(UCSC_ref_addLabel_addType_Chrolist_addCategory.head())

        # UCSC_ref_addLabel_addType_Chrolist_addCategory.to_csv(dest_dir /  f'UCSC_ref_addLabel_addType_chroList_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}_addCategory_usingCDS0{WT_num}D_CDS0{HS_num}D.txt', index=None, sep='\t')
        genes = pd.read_csv(root.parent / 'results' / '21_RNA-seq_with_largestExpression_isoform_add5Category' / f'UCSC_ref_addLabel_addType_chroList_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}_addCategory_usingCDS0{WT_num}D_CDS0{HS_num}D.txt', sep='\t', dtype={'ID': str})
        # print(genes.shape)
        # print(genes.head())
        expressed_raw = genes[genes['expression'] == 'expressed']
        # print(expressed_raw.shape)
        # print(expressed_raw.head())
        expressed = expressed_raw[(expressed_raw['FPKM'] > 0) & (expressed_raw['FPKM_RDS026'] > 0)]
        print(expressed.shape)
        # print(expressed.head())

        expressed_keep = expressed[expressed['category'] == 'keep']
        print(expressed_keep.shape)
        print(expressed_keep.head())
        # print(expressed_keep['gene_symbol'])
        expressed_keep_IDs = map_symbols_to_ids(expressed_keep['gene_symbol'], IDs)
        # print(len(expressed_keep_IDs))
        # print(expressed_keep_IDs)
        expressed_keep_IDs.to_csv(dest_dir / f'_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}_usingCDS0{WT_num}D_CDS0{HS_num}D_expressed_keep.txt',header=None, index=None,sep='\t')
        
        expressed_appear = expressed[expressed['category'] == 'appear']
        # print(expressed_appear.shape)
        # print(expressed_appear.head())
        expressed_appear_IDs = map_symbols_to_ids(expressed_appear['gene_symbol'], IDs)
        # print(len(expressed_appear_IDs))
        # print(expressed_appear_IDs)
        expressed_appear_IDs.to_csv(dest_dir / f'_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}_usingCDS0{WT_num}D_CDS0{HS_num}D_expressed_appear.txt', header=None, index=None, sep='\t')
        
        expressed_disappear = expressed[expressed['category'] == 'disappear']
        # print(expressed_disappear.shape)
        # print(expressed_disappear.head())
        expressed_disappear_IDs = map_symbols_to_ids(expressed_disappear['gene_symbol'], IDs)
        # print(len(expressed_disappear_IDs))
        # print(expressed_disappear_IDs)
        expressed_disappear_IDs.to_csv(dest_dir / f'_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}_usingCDS0{WT_num}D_CDS0{HS_num}D_expressed_disappear.txt',header=None, index=None, sep='\t')

        expressed_merge = expressed[expressed['category'] == 'merge']
        expressed_merge_IDs = map_symbols_to_ids(expressed_merge['gene_symbol'], IDs)
        # print(len(expressed_merge_IDs))
        # print(expressed_merge_IDs)
        expressed_merge_IDs.to_csv(dest_dir / f'_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}_usingCDS0{WT_num}D_CDS0{HS_num}D_expressed_merge.txt',header=None, index=None, sep='\t')


        expressed_split = expressed[expressed['category'] == 'split']
        expressed_split_IDs = map_symbols_to_ids(expressed_split['gene_symbol'], IDs)
        # print(len(expressed_split_IDs))
        # print(expressed_split_IDs)
        expressed_split_IDs.to_csv(dest_dir / f'_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}_usingCDS0{WT_num}D_CDS0{HS_num}D_expressed_split.txt', header=None, index=None, sep='\t')
