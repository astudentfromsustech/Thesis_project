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


def swap_values(row, chrom_order):
    chrom1 = chrom_order[row['chrom1_dm3']]
    chrom2 = chrom_order[row['chrom2_dm3']]

    # Swap the chromosome and position values if chrom1_dm3 > chrom2_dm3
    if chrom1 > chrom2:
        row['chrom1_dm3'], row['chrom2_dm3'] = row['chrom2_dm3'], row['chrom1_dm3']
        row['pos1_dm3'], row['pos2_dm3'] = row['pos2_dm3'], row['pos1_dm3']
    return row

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'pairs'
    anno_dir = root.parent / 'annotations'
    dest_dir = root.parent / 'data' / '0_pairs_from_dm6_to_dm3'
    dest_dir.mkdir(exist_ok=True, parents=True)

    lo = LiftOver(str(anno_dir / f'dm6ToDm3.over.chain.gz'))
    # conditions = ['WT', 'HS']
    conditions = ['WT']
    # pair_nums= [10000000, 20000000, 40000000, 60000000, 80000000, 100000000]
    # pair_nums = [10000000]
    chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
    chrom_order = {'2L': 1, '2R': 2, '3L': 3, '3R': 4, '4': 5, 'X': 6}

    for condition in conditions:
        print(condition)
        # raw_data = pd.read_csv(src_dir / f'{condition}_combined_PNAS.pairs', comment='#', header=None, sep='\s+').head(10)
        raw_data = pd.read_csv(src_dir / f'downsampled_pairs.txt', comment='#', header=None, sep='\s+').head(100)
        raw_data.columns = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type', 'mapq1', 'mapq2']
        raw_data['strand1'] = raw_data['strand1'].replace({'-': '0', '+': '1'})
        raw_data['strand2'] = raw_data['strand2'].replace({'-': '0', '+': '1'})

        print(raw_data.shape)
        print(raw_data.head())
        df = raw_data[['strand1', 'chrom1', 'pos1', 'strand2', 'chrom2', 'pos2', ]]
        df.insert(3, 'frag1', 0)
        df.insert(7, 'frag2', 1)
        print(df.shape)
        print(df.head())

        df['chrom1_dm3'], df['pos1_dm3'] = zip(*df.apply(lambda row: convert_coordinates(row['chrom1'], row['pos1'], lo), axis=1))
        df['chrom2_dm3'], df['pos2_dm3'] = zip(*df.apply(lambda row: convert_coordinates(row['chrom2'], row['pos2'], lo), axis=1))
        # print(df)
        df_dm3 = df[['strand1', 'chrom1_dm3', 'pos1_dm3', 'frag1', 'strand2', 'chrom2_dm3', 'pos2_dm3', 'frag2']]
        # print(df_dm3)
        df_dm3_filterChro_raw = df_dm3[(df_dm3['chrom1_dm3'].isin(chro_list)) & (df_dm3['chrom2_dm3'].isin(chro_list))]
        df_dm3_filterChro = df_dm3_filterChro_raw.copy()
        df_dm3_filterChro['chrom1_dm3'] = df_dm3_filterChro_raw['chrom1_dm3'].str.replace('chr', '')
        df_dm3_filterChro['chrom2_dm3'] = df_dm3_filterChro_raw['chrom2_dm3'].str.replace('chr', '')
        print(df_dm3_filterChro.shape)
        print(df_dm3_filterChro.head())

        df_dm3_filterChro_exchange = df_dm3_filterChro.apply(swap_values, chrom_order=chrom_order, axis=1)

        data_sorted = df_dm3_filterChro_exchange.sort_values(by=['chrom1_dm3', 'chrom2_dm3'])
        data_sorted['pos1_dm3'] = data_sorted['pos1_dm3'].astype(int)
        data_sorted['pos2_dm3'] = data_sorted['pos2_dm3'].astype(int)
        print(data_sorted.shape)
        print(data_sorted.head())
        pair_counts = data_sorted .groupby(['chrom1_dm3', 'chrom2_dm3']).size().reset_index(name='counts')
        print(pair_counts)
        unique_pairs = data_sorted[['chrom1_dm3', 'chrom2_dm3']].drop_duplicates()
        print(unique_pairs)

        # data_sorted.to_csv(dest_dir / f'{condition}_combined_PNAS.pairs.Pre', index=None, header=None, sep=' ')


