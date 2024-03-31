import pathlib as p
import pandas as pd
import numpy as np

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
    src_dir = root.parent / 'data' / '0_pairs_from_dm6_to_dm3_usingChunks'
    anno_dir = root.parent / 'annotations'
    dest_dir = root.parent / 'data' / '1_pairs_from_dm6_to_dm3_usingChunks_combineChunks'
    dest_dir.mkdir(exist_ok=True, parents=True)

    conditions = ['WT', 'HS']
    # conditions = ['HS']

    for condition in conditions:
        print(condition)
        base_filename = f'{condition}_combined_PNAS_chunk*.pairs.Pre'

        data_combined = pd.DataFrame()
        for file_path in src_dir.glob(base_filename):
            print(file_path)
            data = pd.read_csv(file_path, header=None, sep=' ')
            data.columns = ['strand1', 'chrom1_dm3', 'pos1_dm3', 'frag1', 'strand2', 'chrom2_dm3', 'pos2_dm3', 'frag2']
            print(data.shape)
            # print(data.head())
            data_combined = pd.concat([data_combined, data])
            print(data_combined.shape)
        print(data_combined.shape)
        data_combined['pos1_dm3'] = data_combined['pos1_dm3'].astype(int)
        data_combined['pos2_dm3'] = data_combined['pos2_dm3'].astype(int)
        print(data_combined.shape)
        print(data_combined.head())
        data_combined.to_csv(dest_dir / f'{condition}_combined_PNAS_chunksCombined.pairs.Pre', index=None, header=None, sep=' ')


