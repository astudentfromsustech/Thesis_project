import pathlib as p
import pandas as pd

# import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir_cis = root.parent / 'data' / 'combined_cis_cluster_fiterChroms_addLoopSpan'
    dest_dir = root.parent / 'data' / 'combined_cis_cluster_fiterChroms_addLoopSpan_filterPETcount_forIGV'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 5
    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        print(num)
        loops = pd.read_csv(src_dir_cis / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan.txt', sep='\t')

        print(loops.shape)
        print(loops.head())
        loops_fiterCount = loops[loops['count'] >= PETcount_thresh]
        print(loops_fiterCount.shape)
        print(loops_fiterCount.head())
        #
        loops_fiterCount.to_csv(dest_dir / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_forIGV.bedpe', header=None, index=None, sep='\t')





