# import cooler
import pathlib as p
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from matplotlib.colors import LinearSegmentedColormap

from scipy.stats import pearsonr, spearmanr, gaussian_kde, ttest_ind

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)


def call_interval(df, chro, res):
    """
    Adds three columns ('chr', 'start', 'end') to the beginning of the dataframe.
    The 'chr' column is filled with the string 'chr'.
    The 'start' and 'end' columns are filled with numeric ranges starting from 1 to 10000 for the first row,
    increasing by 10000 for each subsequent row for 'start' and by 10000 for 'end', but with a different initial offset.
    """
    # Number of rows in the dataframe
    n = len(df)

    # Generating the 'chr', 'start', and 'end' columns
    chr_column = [chro] * n
    start_column = [1 + res * i for i in range(n)]
    end_column = [res + res * i for i in range(n)]

    # Inserting the new columns at the beginning of the dataframe
    df.insert(0, 'chromosome', chr_column)
    df.insert(1, 'start', start_column)
    df.insert(2, 'end', end_column)

    return df

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent

    chro_list = ['chr2L']
    num_list = ['WT']
    resolution='10k'
    res = 10000
    # chro = 'chrX'
    # num = '01'
    for num in num_list:
        print(num)
        for chro in chro_list:
            print(chro)
            src_dir = root.parent / 'results' / 'results_combined_data' / f'eigen_vector_{chro}_{resolution}_KR_called_by_juicer_tools'
            dest_dir = root.parent / 'results' / 'results_combined_data' / f'0_call_AB_compartment_intervals'
            dest_dir.mkdir(parents=True, exist_ok=True)

            eigens = pd.read_csv(src_dir / f'{num}_combined_eigen_{chro}_{resolution}_KR.txt', header=None, sep=' ').iloc[:, 0].to_frame(name='score')
            print(type(eigens))
            print(eigens.shape)
            print(eigens.head())
            compartment_intervals = call_interval(eigens, chro, res)
            compartment_intervals['score'].fillna(0, inplace=True)
            compartment_intervals['score'] = compartment_intervals['score']*(-1)
            print(compartment_intervals.shape)
            print(compartment_intervals)

            compartment_intervals.to_csv(dest_dir / f'compartment_intervals_{num}_{chro}_res{res}_forJuicebox.bedGraph', header=None, index=None, sep='\t')
