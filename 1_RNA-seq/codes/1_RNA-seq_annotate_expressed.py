import pathlib as p
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def assign_type(df, logFC_thresh, FDR_thresh):
    types = []  # To store the type for each row
    for index, row in df.iterrows():
        if pd.isna(row['logFC']):
            types.append('Non-expressed')
        elif row['logFC'] >= logFC_thresh and row['FDR'] <= FDR_thresh:
            types.append('up_regulated')
        elif ((row['logFC'] > 0) and (row['logFC'] < logFC_thresh)) or ((row['logFC'] > logFC_thresh) and (row['logFC'] > logFC_thresh)):
            types.append('up_NS')
        elif row['logFC'] <= -logFC_thresh  and row['FDR'] < FDR_thresh:
            types.append('dn_regulated')
        else:
            types.append('dn_NS')
    df['type'] = types
    return df

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data'
    anno_dir = root.parent / 'annotations'
    dest_dir = root.parent / 'results'
    dest_dir.mkdir(parents=True, exist_ok=True)

    logFC_thresh = 1
    FDR_thresh = 0.01


    # chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
    # dm3_ref = pd.read_csv(anno_dir / f'BASIC_UCSC_dm3_ref_withHighestIosoform_inWT.txt', sep='\t')[['name','chromosome','strand','start', 'end', 'gene_symbol']]
    # dm3_ref.columns = ['name','chromosome','strand','start', 'end', 'gene_id']
    # print(dm3_ref.shape)
    # print(dm3_ref.head())
    # dm3_filter = dm3_ref[dm3_ref['chromosome'].isin(chro_list)]
    # print(dm3_filter.shape)
    #
    # #
    # genes = pd.read_csv(src_dir / '0_RDS029_RDS026_combined_data.txt', sep='\t')[['gene_id', 'expected_count', 'TPM','FPKM', 'expected_count_RDS026', 'TPM_RDS026', 'FPKM_RDS026']]
    # # print(genes.shape)
    # # print(genes.head())
    # #
    # expressed_genes = pd.read_csv(src_dir / '0_CDS029_CDS026_edgeR_result.txt', sep='\t')[['gene_id', 'logFC', 'logCPM','FDR']]
    # expressed_genes['expressed'] = 'expressed'
    # # print(expressed_genes.shape)
    # # print(expressed_genes.head(10))
    # # #
    # expressed_genes = pd.merge(genes, expressed_genes, on='gene_id', how='left').sort_values(by='FDR', ascending=True, na_position='last').reset_index(drop=True)
    # # print(expressed_genes.shape)
    # # print(expressed_genes.iloc[8200:8280])
    #
    # #
    # expressed_genes_addAnno = pd.merge(expressed_genes, dm3_filter, on='gene_id', how='inner')[
    #     ['gene_id', 'name', 'chromosome', 'strand', 'start', 'end', 'expected_count', 'TPM', 'FPKM', 'expected_count_RDS026', 'TPM_RDS026', 'FPKM_RDS026', 'logFC', 'logCPM', 'FDR', 'expressed']]
    # # print(expressed_genes_addAnno.shape)
    # # print(expressed_genes_addAnno[expressed_genes_addAnno['expressed']=='expressed'].shape)
    # # print(expressed_genes_addAnno.iloc[8000:8100])
    # expressed_genes_addAnno.to_csv(dest_dir / f'1_annotated_all_genes_using_dm3_ref_withHighestIsoform_inWT.txt', index=None, sep='\t')

    genes_anno = pd.read_csv(dest_dir / f'1_annotated_all_genes_using_dm3_ref_withHighestIsoform_inWT.txt', sep='\t')
    # print(genes_anno.shape)
    # print(genes_anno.iloc[8000:8100])
    genes_anno_addType = assign_type(genes_anno, logFC_thresh, FDR_thresh)
    print(genes_anno_addType.shape)
    # print(genes_anno_addType.iloc[1000:2000])
    print(genes_anno_addType['type'].value_counts())

    genes_anno.to_csv(dest_dir / f'1_genes_anno_logFC{logFC_thresh}_FDR{FDR_thresh}.txt', sep='\t', index=None)







