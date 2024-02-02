import pathlib as p
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data'
    anno_dir = root.parent / 'annotations'
    dest_dir = root.parent / 'results'
    dest_dir.mkdir(parents=True, exist_ok=True)

    TPM_thresh = 1

    chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
    dm3_ref = pd.read_csv(anno_dir / f'BASIC_UCSC_dm3_ref_withHighestIosoform_inWT.txt', sep='\t')[['name','chromosome','strand','start', 'end', 'gene_symbol']]
    dm3_ref.columns = ['name','chromosome','strand','start', 'end', 'gene_id']
    # print(dm3_ref.shape)
    # print(dm3_ref.head())
    dm3_filter = dm3_ref[dm3_ref['chromosome'].isin(chro_list)]
    print(dm3_filter.shape)
    print(dm3_filter.head())

    genes = pd.read_csv(src_dir / '0_RDS029_RDS026_combined_data.txt', sep='\t')[['gene_id', 'expected_count', 'TPM','FPKM', 'expected_count_RDS026', 'TPM_RDS026', 'FPKM_RDS026']]
    # print(genes.shape)
    # print(genes.head())

    expressed_genes = pd.read_csv(src_dir / '0_CDS029_CDS026_edgeR_result.txt', sep='\t')[['gene_id', 'logFC', 'logCPM','FDR']]
    # print(expressed_genes.shape)
    # print(expressed_genes.head(10))

    expressed_genes = pd.merge(genes, expressed_genes, on='gene_id', how='left').sort_values(by='FDR', ascending=True, na_position='last')
    print(expressed_genes.shape)
    print(expressed_genes.head())

    expressed_genes_addAnno = pd.merge(expressed_genes, dm3_filter, on='gene_id', how='inner')[
        ['gene_id', 'name', 'chromosome', 'strand', 'start', 'end', 'expected_count', 'TPM', 'FPKM', 'expected_count_RDS026', 'TPM_RDS026', 'FPKM_RDS026', 'logFC', 'logCPM', 'FDR']]
    # print(expressed_genes_addAnno.shape)
    # print(expressed_genes_addAnno.iloc[7900:8100])
    expressed_genes_addAnno.to_csv(dest_dir / f'1_annotated_all_genes_using_dm3_ref_withHighestIsoform_inWT.txt', index=None, sep='\t')






