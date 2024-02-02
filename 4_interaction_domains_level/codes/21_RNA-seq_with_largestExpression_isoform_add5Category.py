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

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data'
    domain_dir = root.parent / 'results' / '14_Feature_of_interactionDomain_Feature4_STEP1_WT_HS_union'
    anno_dir = root.parent / 'annotations'
    dest_dir = root.parent / 'results' / '21_RNA-seq_with_largestExpression_isoform_add5Category'
    dest_dir.mkdir(exist_ok=True, parents=True)

    TPM_thresh = 0
    logFC_thresh = 1
    FDR_thresh = 0.001

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
    # num_list = [['07', '08', '#1087F4', '#99BDCB']]
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
        genes = pd.read_csv(dest_dir / f'UCSC_ref_addLabel_addType_chroList_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}_addCategory_usingCDS0{WT_num}D_CDS0{HS_num}D.txt', sep='\t', dtype={'ID': str})
        # print(genes.shape)
        # print(genes.head())
        expressed_raw = genes[genes['expression'] == 'expressed']
        # print(expressed_raw.shape)
        # print(expressed_raw.head())
        expressed = expressed_raw[(expressed_raw['FPKM'] > 0) & (expressed_raw['FPKM_RDS026'] > 0)]
        print(expressed.shape)
        print(expressed.head())

        expressed_keep = expressed[expressed['category'] == 'keep']
        # print(expressed_keep.shape)
        # print(expressed_keep.head())
        _, p_keep = wilcoxon(expressed_keep['FPKM'], expressed_keep['FPKM_RDS026'])
        # _, p_keep_log = wilcoxon(np.log10(expressed_keep['FPKM']), np.log10(expressed_keep['FPKM_RDS026']))
        # print('p_keep', p_keep)
        # print('p_keep_log', p_keep_log)
        _, p_keep_t = ttest_rel(expressed_keep['FPKM'], expressed_keep['FPKM_RDS026'])
        print('p_keep_t', round(p_keep_t,2))
        
        expressed_appear = expressed[expressed['category'] == 'appear']
        # print(expressed_appear.shape)
        # print(expressed_appear.head())
        _, p_appear = wilcoxon(expressed_appear['FPKM'], expressed_appear['FPKM_RDS026'])
        # print('p_appear', p_appear)
        # _, p_appear_log = wilcoxon(np.log10(expressed_appear['FPKM']), np.log10(expressed_appear['FPKM_RDS026']))
        # print('p_appear_log', p_appear_log)
        _, p_appear_t = ttest_rel(expressed_appear['FPKM'], expressed_appear['FPKM_RDS026'])
        print('p_appear_t', round(p_appear_t, 2))
        
        expressed_disappear = expressed[expressed['category'] == 'disappear']
        # print(expressed_disappear.shape)
        # print(expressed_disappear.head())
        _, p_disappear = wilcoxon(expressed_disappear['FPKM'], expressed_disappear['FPKM_RDS026'])
        # print('p_disappear', p_disappear)
        # _, p_disappear_log = wilcoxon(np.log10(expressed_disappear['FPKM']), np.log10(expressed_disappear['FPKM_RDS026']))
        # print('p_disappear_log', p_disappear_log)
        _, p_disappear_t = ttest_rel(expressed_disappear['FPKM'], expressed_disappear['FPKM_RDS026'])
        print('p_disappear_t', round(p_disappear_t, 2))

        expressed_merge = expressed[expressed['category'] == 'merge']
        # print(expressed_merge.shape)
        # print(expressed_merge.head())
        _, p_merge = wilcoxon(expressed_merge['FPKM'], expressed_merge['FPKM_RDS026'])
        # print('p_merge', p_merge)
        # _, p_merge_log = wilcoxon(np.log10(expressed_merge['FPKM']), np.log10(expressed_merge['FPKM_RDS026']))
        # print('p_merge_log', p_merge_log)
        _, p_merge_t = ttest_rel(expressed_merge['FPKM'], expressed_merge['FPKM_RDS026'])
        print('p_merge_t', round(p_merge_t,2))

        expressed_split = expressed[expressed['category'] == 'split']
        # print(expressed_split.shape)
        # print(expressed_split.head())
        _, p_split = wilcoxon(expressed_split['FPKM'], expressed_split['FPKM_RDS026'])
        # print('p_split', p_split)
        # _, p_split_log = wilcoxon(np.log10(expressed_split['FPKM']), np.log10(expressed_split['FPKM_RDS026']))
        # print('p_split_log', p_split_log)
        _, p_split_t = ttest_rel(expressed_split['FPKM'], expressed_split['FPKM_RDS026'])
        print('p_split_t', round(p_split_t, 2))

        # df = pd.concat([np.log10(expressed_keep['FPKM']), np.log10(expressed_keep['FPKM_RDS026']), np.log10(expressed_appear['FPKM']), np.log10(expressed_appear['FPKM_RDS026']),
        #                 np.log10(expressed_disappear['FPKM']), np.log10(expressed_disappear['FPKM_RDS026']), np.log10(expressed_merge['FPKM']), np.log10(expressed_merge['FPKM_RDS026']),
        #                 np.log10(expressed_split['FPKM']), np.log10(expressed_split['FPKM_RDS026'])], axis=1)
        df = pd.DataFrame({'WT_keep': np.log10(expressed_keep['FPKM']), 'HS_keep': np.log10(expressed_keep['FPKM_RDS026']),
                           'WT_appear': np.log10(expressed_appear['FPKM']), 'HS_appear': np.log10(expressed_appear['FPKM_RDS026']),
                           'WT_disappear': np.log10(expressed_disappear['FPKM']), 'HS_disappear': np.log10(expressed_disappear['FPKM_RDS026']),
                           'WT_merge': np.log10(expressed_merge['FPKM']), 'HS_merge': np.log10(expressed_merge['FPKM_RDS026']),
                           'WT_split': np.log10(expressed_split['FPKM']), 'HS_split': np.log10(expressed_split['FPKM_RDS026'])})
        # print(df.shape)
        # print(df)
        fig, ax = plt.subplots(figsize=(10, 5), dpi=200)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            # spine.set_color('none')
            spine.set_linewidth(1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        colors_pal = [WT_color, HS_color]
        flierprops = dict(marker='o', markerfacecolor='none', markersize=4, markeredgecolor='black',linestyle='none', markeredgewidth=0.6)
        sns.boxplot(data=df, palette=colors_pal, flierprops=flierprops)  # Boxplot
        sns.stripplot(data=df, color='black', size=1, jitter=True)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=15)
        # plt.show()