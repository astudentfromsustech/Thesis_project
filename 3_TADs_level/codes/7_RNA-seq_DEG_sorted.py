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
    # print(annotations.shape)
    # print(annotations.head())

    annotations_brown = pd.read_csv(annotations_dir / f'BASIC_Brown_dm3_ref_longestIsoform.txt', sep='\t')
    # print(annotations_brown.shape)
    # print(annotations_brown.head())


    DEG = pd.read_csv(src_dir / '7_RNA-seq_DEG_for_basicBrowser.txt', sep='\t')
    print(DEG.shape)
    print(DEG.head())

    # DEG_sotred = DEG.sort_values(by='score', key=lambda x: x.abs(), ascending=False)
    # print(DEG_sotred.shape)
    # print(DEG_sotred.tail(20))
    # DEG_sotred.to_csv(src_dir / '7_RNA-seq_DEG_sorted.txt', index=None, sep='\t')

