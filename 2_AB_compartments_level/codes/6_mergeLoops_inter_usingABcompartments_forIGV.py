import pathlib as p
import pandas as pd


import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns
# import typing as t
plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)



if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '6_mergeLoops_inter_usingABcompartments'
    dest_dir = root.parent / 'results' / '6_mergeLoops_inter_usingABcompartments_forIGV'
    dest_dir.mkdir(parents=True, exist_ok=True)

    PETcount_thresh = 1

    merged_loop_countThresh = 2
    # merged_loop_span_thresh = 500000


    inter_type = ['inter_B', 'inter_A']

    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        print(num)
        loops = pd.read_csv(src_dir / f'CDS0{num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM.bedpe', sep='\t')
        print(loops.shape)
        print(loops.head())
        loops.to_csv(dest_dir / f'CDS0{num}D_inter_merged_loops_merged_loop_countThresh{merged_loop_countThresh}_addDomainDistance_addRPM_forIGV.bedpe',header=None, index=None, sep='\t')
















