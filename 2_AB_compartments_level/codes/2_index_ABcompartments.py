import cooler
import pathlib as p
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from hmmlearn import hmm

np.set_printoptions(threshold=np.inf, linewidth=np.inf)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

def find_intersected_intervals(df1, df2):
    intersected_intervals = []
    # Group by chromosome
    for chromosome in df1['chromosome'].unique():
        # Filter intervals by chromosome
        intervals_df1 = df1[df1['chromosome'] == chromosome]
        intervals_df2 = df2[df2['chromosome'] == chromosome]
        # Find intersecting intervals
        for _, row1 in intervals_df1.iterrows():
            for _, row2 in intervals_df2.iterrows():
                # Check if intervals intersect
                if row1['start'] <= row2['end'] and row2['start'] <= row1['end']:
                    # Calculate the intersection
                    start = max(row1['start'], row2['start'])
                    end = min(row1['end'], row2['end'])
                    intersected_intervals.append({
                        'chromosome': chromosome,
                        'start': start,
                        'end': end
                    })
    # Convert the list of dictionaries to a DataFrame
    intersected_df = pd.DataFrame(intersected_intervals)
    return intersected_df

def add_aid_column(df):
    # Generate the 'A_ID' values starting from A0001
    a_id_values = ['A' + str(i).zfill(4) for i in range(1, len(df) + 1)]
    # Add the 'A_ID' column to the dataframe
    df['A_ID'] = a_id_values
    return df

def add_bid_column(df):
    # Generate the 'A_ID' values starting from A0001
    b_id_values = ['B' + str(i).zfill(4) for i in range(1, len(df) + 1)]
    # Add the 'A_ID' column to the dataframe
    df['B_ID'] = b_id_values
    return df

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent

    # bin_size = 5000
    # dir_num = '5k'

    bin_size = 10000
    dir_num = '10k'

    src_dir = root.parent / 'results' / f'1_ABcompartments_classification'
    dest_dir = root.parent / 'results' / '2_index_ABcompartments'
    dest_dir.mkdir(parents=True, exist_ok=True)

    WT_A = pd.read_csv(src_dir / f'WT_A_merged_resolution{bin_size}.bed', header=None, sep='\t')
    WT_A.columns = ['chromosome', 'start', 'end']
    print(WT_A.shape)
    print(WT_A.head(10))

    HS_A = pd.read_csv(src_dir / f'HS_A_merged_resolution{bin_size}.bed', header=None, sep='\t')
    HS_A.columns = ['chromosome', 'start', 'end']
    print(HS_A.shape)
    print(HS_A.head(10))

    intersected_A = find_intersected_intervals(WT_A, HS_A)
    print(intersected_A.shape)
    print(intersected_A.head())

    intersected_A_addID = add_aid_column(intersected_A)
    print(intersected_A_addID.shape)
    print(intersected_A_addID)
    intersected_A_addID.to_csv(dest_dir / f'intersected_A_addID_resolution{bin_size}.bed', header=None, index=None, sep='\t')

    WT_B = pd.read_csv(src_dir / f'WT_B_merged_resolution{bin_size}.bed', header=None, sep='\t')
    WT_B.columns = ['chromosome', 'start', 'end']
    print(WT_B.shape)
    print(WT_B.head(10))
    HS_B = pd.read_csv(src_dir / f'HS_B_merged_resolution{bin_size}.bed', header=None, sep='\t')
    HS_B.columns = ['chromosome', 'start', 'end']
    print(HS_B.shape)
    print(HS_B.head(10))
    intersected_B = find_intersected_intervals(WT_B, HS_B)
    intersected_B_addID = add_bid_column(intersected_B)
    print(intersected_B_addID.shape)
    print(intersected_B_addID)
    intersected_B_addID.to_csv(dest_dir / f'intersected_B_addID_resolution{bin_size}.bed', header=None, index=None, sep='\t')



