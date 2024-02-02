import pathlib as p
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')

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

def promoter_position(row):
    if row['strand'] == '+':
        return row['start'] - 2000, row['start']
    else:
        return row['end'], row['end'] + 2000

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data'
    anno_dir = root.parent / 'annotations'

    TPM_thresh = 0
    logFC_thresh = 1
    FDR_thresh = 0.001



    ## get expressed genes
    expressed_edgeR = pd.read_csv(src_dir/ 'RNA-seq' / 'demo8_expressed_data_RSEM_edegR.txt', sep='\t')[['gene_id', 'TPM', 'TPM_RDS026', 'FPKM', 'FPKM_RDS026','logFC', 'FDR']]
    print(expressed_edgeR.shape)
    print(expressed_edgeR.head())
    expressed_edgeR_filterTPM = expressed_edgeR[(expressed_edgeR['TPM'] >= TPM_thresh) | (expressed_edgeR['TPM_RDS026'] >= TPM_thresh)]
    print(expressed_edgeR_filterTPM.shape)
    print(expressed_edgeR_filterTPM.head())
    expressed_edgeR_filterTPM_filterHis = expressed_edgeR_filterTPM[~expressed_edgeR_filterTPM['gene_id'].str.contains(r'^His\w+:CG\d+', na=False)]
    print(expressed_edgeR_filterTPM_filterHis.shape)
    print(expressed_edgeR_filterTPM_filterHis.head())
    #
    UCSC_ref = pd.read_csv(anno_dir / '1_UCSC_ref_withHighestIosoform_inWT.txt', sep='\t').drop(columns=['mRNA', 'exons', 'CDS_start', 'CDS_end'])
    # UCSC_ref = pd.read_csv(anno_dir / 'BASIC_UCSC_dm3_ref_sortedBymRNA_longestIsoform.txt', sep='\t').drop(columns=['mRNA', 'exons', 'CDS_start', 'CDS_end'])
    UCSC_ref.sort_values(by=['chromosome', 'start'], inplace=True)
    print(UCSC_ref.shape)
    print(UCSC_ref.head())
    UCSC_ref_addLabel = UCSC_ref.copy()
    UCSC_ref_addLabel[['P_left', 'P_right']] = UCSC_ref.apply(lambda row: promoter_position(row), axis=1, result_type='expand')
    UCSC_ref_addLabel['expression'] = UCSC_ref['gene_symbol'].apply(lambda x: 'expressed' if x in expressed_edgeR_filterTPM_filterHis['gene_id'].values else 'unexpressed')
    print(UCSC_ref_addLabel.shape)
    print(UCSC_ref_addLabel.head())
    UCSC_ref_addLabel[['tss', 'tes']] = UCSC_ref.apply(lambda row: (row['start'], row['end']) if row['strand'] == '+' else (row['end'], row['start']), axis=1, result_type='expand')
    UCSC_ref_addLabel[['TPM', 'TPM_RDS026', 'FPKM', 'FPKM_RDS026', 'logFC', 'FDR']] = UCSC_ref_addLabel.apply(lambda row: get_expression_data(row, expressed_edgeR_filterTPM_filterHis), axis=1, result_type='expand')
    print(UCSC_ref_addLabel.shape)
    print(UCSC_ref_addLabel.head())
    #
    UCSC_ref_addLabel_addType = UCSC_ref_addLabel.copy()
    UCSC_ref_addLabel_addType.loc[(UCSC_ref_addLabel_addType['logFC'] > logFC_thresh) & (UCSC_ref_addLabel_addType['FDR'] < FDR_thresh), 'type'] = 'up_regulated'
    UCSC_ref_addLabel_addType.loc[ (UCSC_ref_addLabel_addType['logFC'] < -logFC_thresh) & (UCSC_ref_addLabel_addType['FDR'] < FDR_thresh), 'type'] = 'down_regulated'
    UCSC_ref_addLabel_addType.loc[((UCSC_ref_addLabel_addType['logFC'] < logFC_thresh) & (UCSC_ref_addLabel_addType['logFC'] > 0)) | ((UCSC_ref_addLabel_addType['logFC'] > logFC_thresh) & (UCSC_ref_addLabel_addType['FDR'] > FDR_thresh)), 'type'] = 'pos_NS'
    UCSC_ref_addLabel_addType.loc[((UCSC_ref_addLabel_addType['logFC'] < 0) & (UCSC_ref_addLabel_addType['logFC'] > -logFC_thresh)) | ( (UCSC_ref_addLabel_addType['logFC'] < -logFC_thresh) & (UCSC_ref_addLabel_addType['FDR'] > FDR_thresh)), 'type'] = 'neg_NS'
    chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
    UCSC_ref_addLabel_addType_ChroList = UCSC_ref_addLabel_addType[UCSC_ref_addLabel_addType['chromosome'].isin(chro_list)]
    print(UCSC_ref_addLabel_addType_ChroList.shape)
    print(UCSC_ref_addLabel_addType_ChroList.head())
    UCSC_ref_addLabel_addType_ChroList.to_csv(src_dir/ 'RNA-seq' / f'UCSC_ref_addLabel_addType_chroList_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}.txt', index=None, sep='\t')
    # UCSC_ref_addLabel_addType = pd.read_csv(src_dir / 'RNA-seq' / f'UCSC_ref_addLabel_addType_TPM{TPM_thresh}_FDR{FDR_thresh}_logFC{logFC_thresh}.txt', sep='\t')
    # print(UCSC_ref_addLabel_addType.shape)
    # print(UCSC_ref_addLabel_addType.head(100))
    # expressed = UCSC_ref_addLabel_addType[UCSC_ref_addLabel_addType['expression'] == 'expressed']
    # print(expressed.shape)
    # print(expressed.head())
    # up_DEGs = UCSC_ref_addLabel_addType[UCSC_ref_addLabel_addType['type'] == 'up_regulated']
    # print(up_DEGs.shape)
    # print(up_DEGs.head())
    # down_DEGs = UCSC_ref_addLabel_addType[UCSC_ref_addLabel_addType['type'] == 'down_regulated']
    # print(down_DEGs.shape)
    # print(down_DEGs.head())