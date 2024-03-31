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
    dest_dir = root.parent / 'data' / '0_pairs_from_dm6_to_dm3_usingChunks'
    dest_dir.mkdir(exist_ok=True, parents=True)

    # Initialize LiftOver object
    lo = LiftOver(str(anno_dir / 'dm6ToDm3.over.chain.gz'))
    conditions = ['WT', 'HS']
    chro_list = ['2L', '2R', '3L', '3R', '4', 'X']

    chunk_size = 1000000  # Define a suitable chunk size

    for condition in conditions:
        print(f'Processing condition: {condition}')
        filepath = src_dir / f'{condition}_combined_PNAS.pairs'

        # Initialize a counter for chunk naming
        chunk_counter = 1

        # Process each chunk
        for chunk_raw in pd.read_csv(filepath, comment='#', header=None, sep='\s+', chunksize=chunk_size):
            chunk_raw.columns = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type', 'mapq1','mapq2']
            chunk_raw['strand1'] = chunk_raw['strand1'].replace({'-': '0', '+': '1'})
            chunk_raw['strand2'] = chunk_raw['strand2'].replace({'-': '0', '+': '1'})

            chunk_tmp = chunk_raw[['strand1', 'chrom1', 'pos1', 'strand2', 'chrom2', 'pos2']]
            chunk = chunk_tmp.copy()
            chunk[['chrom1', 'pos1']] = chunk_tmp.apply(lambda row: convert_coordinates(row['chrom1'], row['pos1'], lo),
                                                    axis=1, result_type='expand')
            chunk[['chrom2', 'pos2']] = chunk_tmp.apply(lambda row: convert_coordinates(row['chrom2'], row['pos2'], lo),
                                                    axis=1, result_type='expand')

            # Filter based on chromosome list and adjust chromosome names
            chunk_filtered_raw = chunk[(chunk['chrom1'].str[3:].isin(chro_list)) & (chunk['chrom2'].str[3:].isin(chro_list))]
            chunk_filtered = chunk_filtered_raw.copy()
            chunk_filtered['chrom1'] = chunk_filtered_raw['chrom1'].str[3:]
            chunk_filtered['chrom2'] = chunk_filtered_raw['chrom2'].str[3:]
            chunk_filtered_sorted = chunk_filtered[chunk_filtered['chrom1'] <= chunk_filtered['chrom2']]
            chunk_filtered_sorted.insert(3, 'frag1', 0)
            chunk_filtered_sorted.insert(7, 'frag2', 1)

            print(chunk_filtered_sorted.shape)
            print(chunk_filtered_sorted.head())
            # Define the filename for the processed chunk
            chunk_filename = f'{condition}_combined_PNAS_chunk{chunk_counter}.pairs.Pre'
            chunk_path = dest_dir / chunk_filename

            # Save the processed chunk
            chunk_filtered_sorted.to_csv(chunk_path, index=False, header=None, sep=' ')
            print(f'Chunk {chunk_counter} processed and saved to {chunk_path}')
            chunk_counter += 1

        print(f'All chunks for condition {condition} processed.')


