import pathlib as p
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def determine_label(row, intervals):
    tss = row['TSS']
    for _, interval in intervals.iterrows():
        if interval['Start'] <= tss <= interval['End']:
            return interval['Label']

def calculate_promoter(row, promoter_length=1000):
    if row['strand'] == '+':
        return row['start'] - promoter_length, row['start']
    else:
        return row['end'], row['end'] + promoter_length

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'RNA-seq'
    annotations_dir =root.parent / 'annotations'

    annotations = pd.read_csv(annotations_dir / f'BASIC_UCSC_dm3_ref_sortedBymRNA_longestIsoform.txt', sep='\t')
    # print(annotations.shape)
    # print(annotations.head())

    annotations_brown = pd.read_csv(annotations_dir / f'BASIC_Brown_dm3_ref_longestIsoform.txt', sep='\t')
    # print(annotations_brown.shape)
    # print(annotations_brown.head())


    expressed = pd.read_csv(src_dir / 'demo8_expressed_data_RSEM_edegR.txt', sep='\t')
    # print(expressed.shape)
    # print(expressed.head())

    expressed_addLabels = expressed.copy()
    expressed_addLabels.loc[(expressed_addLabels['logFC'] > 1) & (expressed_addLabels['FDR'] < 0.001), 'type'] = 'up_regulated'
    expressed_addLabels.loc[(expressed_addLabels['logFC'] < -1) & (expressed_addLabels['FDR'] < 0.001), 'type'] = 'down_regulated'
    expressed_addLabels.loc[((expressed_addLabels['logFC'] < 1) & (expressed_addLabels['logFC'] > 0)) | ((expressed_addLabels['logFC'] > 1) & (expressed_addLabels['FDR'] > 0.001)), 'type'] = 'pos_NS'
    expressed_addLabels.loc[((expressed_addLabels['logFC'] < 0) & (expressed_addLabels['logFC'] > -1)) | ((expressed_addLabels['logFC'] < -1) & (expressed_addLabels['FDR'] > 0.001)) , 'type'] = 'neg_NS'
    print(expressed_addLabels.shape)
    print(expressed_addLabels.head())
    print(expressed_addLabels['type'].value_counts())

    result = pd.merge(expressed_addLabels, annotations[['gene_symbol', 'chromosome', 'strand', 'start', 'end']],
                 left_on='gene_id', right_on='gene_symbol', how='left')

    # Drop the 'gene_symbol' column from the result as it's redundant
    result.drop('gene_symbol', axis=1, inplace=True)
    print(result.shape)
    print(result.head())

    #
    # result_filterNa = result[result['start'].isna()]
    # print(result_filterNa.shape)
    # print(result_filterNa)
    #
    result = pd.merge(result, annotations_brown[['Gene Symbol', 'Chromosome', 'Strand','Start', 'End']],
                 left_on='gene_id', right_on='Gene Symbol', how='left')
    print(result.shape)
    print(result.head())
    #
    na_rows = result['chromosome'].isna()
    result.loc[na_rows, 'chromosome'] = result.loc[na_rows, 'Chromosome']
    result.loc[na_rows, 'strand'] = result.loc[na_rows, 'Strand']
    result.loc[na_rows, 'start'] = result.loc[na_rows, 'Start']
    result.loc[na_rows, 'end'] = result.loc[na_rows, 'End']
    print(result.shape)
    print(result.head())
    #
    # print(result.loc[na_rows].shape)
    # print(result.loc[na_rows])

    # result_filterNa = result[result['start'].isna()]
    # print(result_filterNa.shape)
    # print(result_filterNa)
    #
    result.dropna(subset=['chromosome'], inplace=True)
    print(result.shape)
    print(result.head())
    # #
    # result_filterNa = result[result['start'].isna()]
    # print(result_filterNa.shape)
    # print(result_filterNa)

    result_arrange = result[['gene_id', 'chromosome', 'strand', 'start', 'end', 'effective_length',
                             'expected_count', 'TPM', 'FPKM', 'expected_count_RDS026',
                             'TPM_RDS026', 'FPKM_RDS026', 'logFC', 'logCPM', 'PValue', 'FDR', 'type']]
    print(result_arrange.shape)
    print(result_arrange.head())


    result_arrange_addPromoter = result_arrange.copy()
    result_arrange_addPromoter[['promoter_start', 'promoter_end']] = result_arrange.apply(
        lambda row: pd.Series(calculate_promoter(row)), axis=1)
    print(result_arrange_addPromoter.shape)
    print(result_arrange_addPromoter.head())

    result_arrange_addPromoter_sortedChroms = result_arrange_addPromoter.sort_values(by=['chromosome', 'start'])
    print(result_arrange_addPromoter_sortedChroms.shape)
    print(result_arrange_addPromoter_sortedChroms.head(10))

    result_arrange_addPromoter_sortedChroms.to_csv(src_dir / '8_step1_expressed_genes_allInfo_includingPromoter_sortedByChrom.txt', index=None, sep='\t')
