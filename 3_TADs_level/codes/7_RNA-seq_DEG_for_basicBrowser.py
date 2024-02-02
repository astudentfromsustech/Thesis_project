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

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'RNA-seq'
    annotations_dir =root.parent / 'annotations'

    annotations = pd.read_csv(annotations_dir / f'BASIC_UCSC_dm3_ref_sortedBymRNA_longestIsoform.txt', sep='\t')
    print(annotations.shape)
    print(annotations.head())

    annotations_brown = pd.read_csv(annotations_dir / f'BASIC_Brown_dm3_ref_longestIsoform.txt', sep='\t')
    print(annotations_brown.shape)
    print(annotations_brown.head())





    expressed = pd.read_csv(src_dir / 'demo8_expressed_data_RSEM_edegR.txt', sep='\t')
    print(expressed.shape)
    print(expressed.head())

    expressed_addLabels = expressed.copy()
    expressed_addLabels.loc[(expressed_addLabels['logFC'] > 1) & (expressed_addLabels['FDR'] < 0.001), 'type'] = 'up_regulated'
    expressed_addLabels.loc[(expressed_addLabels['logFC'] < -1) & (expressed_addLabels['FDR'] < 0.001), 'type'] = 'down_regulated'
    expressed_addLabels.loc[((expressed_addLabels['logFC'] < 1) & (expressed_addLabels['logFC'] > 0)) | ((expressed_addLabels['logFC'] > 1) & (expressed_addLabels['FDR'] > 0.001)), 'type'] = 'pos_NS'
    expressed_addLabels.loc[((expressed_addLabels['logFC'] < 0) & (expressed_addLabels['logFC'] > -1)) | ((expressed_addLabels['logFC'] < -1) & (expressed_addLabels['FDR'] > 0.001)) , 'type'] = 'neg_NS'
    print(expressed_addLabels.shape)
    print(expressed_addLabels.head())
    print(expressed_addLabels['type'].value_counts())

    result = pd.merge(expressed_addLabels, annotations[['gene_symbol', 'chromosome', 'start', 'end']],
                 left_on='gene_id', right_on='gene_symbol', how='left')

    # Drop the 'gene_symbol' column from the result as it's redundant
    result.drop('gene_symbol', axis=1, inplace=True)
    print(result.shape)
    print(result.head())


    conditions = [
        (expressed_addLabels['type'] == 'up_regulated'),
        (expressed_addLabels['type'] == 'down_regulated'),
        (expressed_addLabels['type'] == 'pos_NS'),
        (expressed_addLabels['type'] == 'neg_NS')
    ]
    colors = [
        '255,0,0',  # Red for up_regulated
        '0,0,255',  # Blue for down_regulated
        '255,182,193',  # Light Red for pos_NS
        '173,216,230'  # Light Blue for neg_NS
    ]
    result['color'] =  np.select(conditions, colors, default='unknown')
    print(result.shape)
    print(result.head())
    # print(result['color'].value_counts())


    # result_filterNa = result[result['start'].isna()]
    # print(result_filterNa.shape)
    # print(result_filterNa)

    result = pd.merge(result, annotations_brown[['Gene Symbol', 'Chromosome', 'Start', 'End']],
                 left_on='gene_id', right_on='Gene Symbol', how='left')
    print(result.shape)
    print(result.head())

    na_rows = result['chromosome'].isna()
    result.loc[na_rows, 'chromosome'] = result.loc[na_rows, 'Chromosome']
    result.loc[na_rows, 'start'] = result.loc[na_rows, 'Start']
    result.loc[na_rows, 'end'] = result.loc[na_rows, 'End']
    print(result.shape)
    print(result.head())

    result.dropna(subset=['chromosome'], inplace=True)
    print(result.shape)
    print(result.head())

    # result_filterNa = result[result['start'].isna()]
    # print(result_filterNa.shape)
    # print(result_filterNa)

    genes_Browser = pd.DataFrame(np.zeros((result.shape[0], 7)))
    genes_Browser.columns = ['#Chromosome', 'Start', 'End', 'type', 'score', 'name', 'color']

    genes_Browser['#Chromosome'] = result['chromosome']
    genes_Browser['Start'] = result['start'].astype(int)
    genes_Browser['End'] = result['end'].astype(int)
    genes_Browser['type'] = result['type']
    genes_Browser['score'] = result['logFC'].round(2)
    genes_Browser['color'] = result['color']
    genes_Browser['name'] = result['gene_id']
    print(genes_Browser.shape)
    print(genes_Browser.head())

    genes_Browser.to_csv(src_dir / '7_RNA-seq_DEG_for_basicBrowser.txt', index=None, sep='\t')
