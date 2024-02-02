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
    src_dir_cis = root.parent / 'data' / 'combined_cis_cluster'
    dest_dir = root.parent / 'data' / 'combined_cis_cluster_fiterChroms_addLoopSpan'
    dest_dir.mkdir(parents=True, exist_ok=True)

    chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']

    num_list = ['01', '02', '04', '05', '07', '08']
    for num in num_list:
        print(num)
        loops = pd.read_csv(src_dir_cis / f'CDS0{num}D_combined_clusters_cis.txt', header=None, sep='\t')
        loops.columns = ['chromosome1', 'start1', 'end1', 'chromosome2', 'start2', 'end2', 'count']
        print(loops.shape)
        print(loops.head())
        loops_fiterChroms = loops[loops['chromosome1'].isin(chro_list)]
        print(loops_fiterChroms.shape)
        loops_fiterChroms_addLoopSpan = loops_fiterChroms.copy()
        loops_fiterChroms_addLoopSpan['loop_span'] = loops_fiterChroms.apply(lambda x: int(0.5*(x['end2']+x['start2']-x['end1']-x['start1'])), axis=1)
        print(loops_fiterChroms_addLoopSpan.shape)
        print(loops_fiterChroms_addLoopSpan.head())

        loops_fiterChroms_addLoopSpan.to_csv(dest_dir / f'CDS0{num}D_combined_clusters_cis_filterChrom_addLoopSpan.txt', index=None, sep='\t')





