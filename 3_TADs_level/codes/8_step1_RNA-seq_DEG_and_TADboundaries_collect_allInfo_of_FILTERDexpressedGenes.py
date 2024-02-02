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

# from scipy.stats import pearsonr, spearmanr, gaussian_kde

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'RNA-seq'
    genes = pd.read_csv(src_dir / '8_step1_expressed_genes_allInfo_includingPromoter_sortedByChrom_promoterLen_minus1k.txt', sep='\t')
    # print(genes.shape)
    # print(genes.head())
    genes_filter = genes[(genes['FPKM'] > 5) | (genes['FPKM_RDS026'] > 5)]
    print(genes_filter.shape)
    print(genes_filter.iloc[100:200])

    # genes = genes[genes['FPKM'] > 0]
    # print(genes.shape)
    # print(genes.head())
    #
    # hist, bins = np.histogram(np.log10(genes['FPKM']), bins=100, density=True)
    # cumulative = np.cumsum(hist) / np.sum(hist)
    #
    # fig, ax = plt.subplots(figsize=(5, 5))
    # ax.plot(bins[:-1], cumulative)
    #
    # # Setting labels and title
    # ax.set_title('Cumulative Density Plot')
    # ax.set_xlabel('Value')
    # ax.set_ylabel('Cumulative Density')
    #
    # # Display the plot
    # plt.show()