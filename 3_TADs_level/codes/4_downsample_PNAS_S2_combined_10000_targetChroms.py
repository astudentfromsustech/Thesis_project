import cooler
import pathlib as p
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from matplotlib.colors import LinearSegmentedColormap

# from scipy.stats import pearsonr, spearmanr, gaussian_kde

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def get_total_read_pairs(c):
    # c = cooler.Cooler(cool_file)

    # Fetch the entire matrix and sum all the values
    total_reads = sum(c.matrix(balance=False).fetch(chrom).sum().sum() for chrom in c.chromnames)
    return total_reads

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'PNAS_data_S2_downsample' / 'S2_cool_10k_targetChroms'

    contact_num = 10000000
    sample = 'HS'
    dest_dir = root.parent / 'data' / 'PNAS_data_S2_downsample' / f'S2_cool_10k_targetChroms_downsample_{contact_num}'
    dest_dir.mkdir(parents=True, exist_ok=True)


    c = cooler.Cooler(str(src_dir / f'{sample}_combined_10000_targetChroms.cool'))
    fraction = contact_num / get_total_read_pairs(c)
    bins = c.bins()[:]
    pixels_df = c.pixels()[:]
    pixels_df['count'] = (pixels_df['count'] * fraction).astype(int)


    cooler.create_cooler(cool_uri=str(dest_dir / f'{sample}_combined_10k_targetChroms_downsample_{contact_num}.cool'), bins=bins, pixels=pixels_df)
