import cooler
import pathlib as p
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')

from scipy.stats import pearsonr, spearmanr, gaussian_kde

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)


def filter_zero_indices(list1, list2):
    """
    Filters out indices where either of the lists has a zero.

    Args:
    - list1 (list): First input list.
    - list2 (list): Second input list, should be of the same length as list1.

    Returns:
    - tuple of two lists: Filtered versions of list1 and list2.
    """

    # Check if lists are of the same length
    if len(list1) != len(list2):
        raise ValueError("The input lists must have the same length.")

    # Get indices where neither list has a zero
    non_zero_indices = [i for i, (x, y) in enumerate(zip(list1, list2)) if x != 0 and y != 0]

    # Use list comprehension to filter the lists based on the non-zero indices
    list1_filtered = [list1[i] for i in non_zero_indices]
    list2_filtered = [list2[i] for i in non_zero_indices]

    return list1_filtered, list2_filtered

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'cool_10k_targetChroms_normalizeSmallest'
    dest_dir = root.parent / 'results' / 'result_0_scatterPlot_of_twoReplicates_allChroms_noLabels_noTickLabels_withColorCode'
    dest_dir.mkdir(parents=True, exist_ok=True)
    # print(src_dir)

    groups = [['01', '10', 'H3K27me3_WT', '#4D7731'], ['02', '11', 'H3K27me3_HS', '#98A246'], ['04', '13', 'H3K27ac_WT', '#F50CB8'],
              ['05', '14', 'H3K27ac_HS', '#743CA6'], ['07', '16', 'RNAPII_WT', '#1087F4'], ['08', '17', 'RNAPII_HS', '#99BDCB']]

    for group in groups:
        rep1 = group[0]
        rep2 = group[1]
        group_name = group[2]
        color_code = group[3]
        print(group_name)
        print(rep1)
        print(rep2)
        c1 = cooler.Cooler(str(src_dir / f'CDS0{rep1}D_10000_targetChroms_normalizeSmallest.cool'))
        c2 = cooler.Cooler(str(src_dir / f'CDS0{rep2}D_10000_targetChroms_normalizeSmallest.cool'))

        matrix1_flat_raw = c1.matrix(balance=False)[:].flatten()
        matrix2_flat_raw = c2.matrix(balance=False)[:].flatten()
        # print(matrix1_flat_raw.shape)
        # print(matrix2_flat_raw.shape)
        #
        #
        matrix1_flat, matrix2_flat = filter_zero_indices(matrix1_flat_raw, matrix2_flat_raw)
        # print(len(matrix1_flat))
        # print(len(matrix2_flat))
        corr_pearson, p_pearson = pearsonr(matrix1_flat, matrix2_flat)
        # corr_spearman, p_spearman = spearmanr(matrix1_flat, matrix2_flat)

        fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(1)

        plt.scatter(np.log10(matrix1_flat), np.log10(matrix2_flat), alpha=0.3, s=1.5, c=color_code)
        # plt.scatter(np.log2(matrix1_flat), np.log2(matrix2_flat), alpha=0.3, s=1.5, c='black')
        # sc = ax.scatter(x, y, c=z, s=1.5, alpha=0.3, cmap='inferno')

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels('')
        ax.set_yticks([1, 2, 3])
        ax.set_yticklabels('')

        # ax.set_xlabel(f"{group_name}_rep1", fontsize=10)
        # ax.set_ylabel(f"{group_name}_rep2", fontsize=10)
        ax.set_xlabel(f"")
        ax.set_ylabel(f"")

        ax.set_title(f"")

        plt.savefig(dest_dir/f'{group_name}_correlation_twoReps_allChroms_pearsonCorr{corr_pearson:.2f}_noLabels_withColorCode.png', dpi=300)
        # plt.show()

