import pathlib as p
import pandas as pd
import numpy as np

from pyliftover import LiftOver

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def convert_coordinates(chrom, pos, lo):
    converted = lo.convert_coordinate(chrom, pos-1)  # Convert from 1-based to 0-based
    if converted:
        # Convert back to 1-based position and return
        return converted[0][0], converted[0][1] + 1
    else:
        # Return NaN for unmapped coordinates
        return np.nan, np.nan


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'pairs'
    anno_dir = root.parent / 'annotations'
    dest_dir = root.parent / 'data' / '0_pairs_from_dm6_to_dm3'
    dest_dir.mkdir(exist_ok=True, parents=True)

    lo = LiftOver(str(anno_dir / f'dm6ToDm3.over.chain.gz'))
    conditions = ['WT', 'HS']
    # conditions = ['WT']
    chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']


    for condition in conditions:
        print(condition)
        raw_data = pd.read_csv(src_dir / f'{condition}_combined_PNAS.pairs', comment='#', header=None, sep='\s+')
        # raw_data = pd.read_csv(src_dir / f'downsampled_pairs.txt', comment='#', header=None, sep='\s+').head(10)
        raw_data.columns = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type', 'mapq1', 'mapq2']
        raw_data['strand1'] = raw_data['strand1'].replace({'-': '0', '+': '1'})
        raw_data['strand2'] = raw_data['strand2'].replace({'-': '0', '+': '1'})

        print(raw_data.shape)
        print(raw_data.head())
        df_raw = raw_data[['strand1', 'chrom1', 'pos1', 'strand2', 'chrom2', 'pos2']]

        df_dm3 = df_raw.copy()
        df_dm3['chrom1'], df_dm3['pos1'] = zip(*df_raw.apply(lambda row: convert_coordinates(row['chrom1'], row['pos1'], lo), axis=1))
        df_dm3['chrom2'], df_dm3['pos2'] = zip(*df_raw.apply(lambda row: convert_coordinates(row['chrom2'], row['pos2'], lo), axis=1))
        print(df_dm3.shape)
        print(df_dm3.head())
        # df_dm3 = df[['strand1', 'chrom1_dm3', 'pos1_dm3', 'frag1', 'strand2', 'chrom2_dm3', 'pos2_dm3', 'frag2']]
        # # print(df_dm3)
        df_dm3_filterChro_raw = df_dm3[(df_dm3['chrom1'].isin(chro_list)) & (df_dm3['chrom2'].isin(chro_list))]
        df_dm3_filterChro = df_dm3_filterChro_raw.copy()
        df_dm3_filterChro['chrom1'] = df_dm3_filterChro_raw['chrom1'].str.replace('chr', '')
        df_dm3_filterChro['chrom2'] = df_dm3_filterChro_raw['chrom2'].str.replace('chr', '')
        print(df_dm3_filterChro.shape)
        print(df_dm3_filterChro.head())
        #
        df_dm3_filterChroOrder = df_dm3_filterChro[df_dm3_filterChro['chrom1'] <= df_dm3_filterChro['chrom2']]

        df_dm3_filterChroOrder.insert(3, 'frag1', 0)
        df_dm3_filterChroOrder.insert(7, 'frag2', 1)
        print(df_dm3_filterChroOrder.shape)
        print(df_dm3_filterChroOrder.head())

        df_dm3_filterChroOrder.to_csv(dest_dir / f'{condition}_combined_PNAS.pairs.Pre', index=None, header=None, sep=' ')


