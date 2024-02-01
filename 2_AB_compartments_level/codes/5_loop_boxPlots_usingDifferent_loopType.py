import pathlib as p
import pandas as pd

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)



if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '4_annotateLoops'
    dest_dir = root.parent / 'results' / '5_loop_boxPlots_usingDifferent_loopType'
    dest_dir.mkdir(parents=True, exist_ok=True)

    bin_size = 10000

    PETcount_thresh = 1
    # loop_span_thresh = 500000

    scatterPlot_PETthresh = 2



    # num_list = [['01', '02'], ['04', '05'], ['07', '08']]
    num_list = [['04', '05']]
    for num_pair in num_list:
        WT_num = num_pair[0]
        HS_num = num_pair[1]
        print(num_pair)
        print(WT_num)
        print(HS_num)

        WT_loops_anno_addType_raw = pd.read_csv(src_dir / f'ANNO_CDS0{WT_num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone_addType.txt', sep='\t')
        WT_loops_anno_addType = WT_loops_anno_addType_raw[WT_loops_anno_addType_raw['count'] >= scatterPlot_PETthresh]
        print(WT_loops_anno_addType.shape)
        WT_intra_A = WT_loops_anno_addType[WT_loops_anno_addType['type'] == 'intra_A']
        print(WT_intra_A.shape)
        print(WT_intra_A.head())

        WT_inter_A = WT_loops_anno_addType[WT_loops_anno_addType['type'] == 'inter_A']
        print(WT_inter_A.shape)
        print(WT_inter_A.head())

        WT_intra_B = WT_loops_anno_addType[WT_loops_anno_addType['type'] == 'intra_B']
        print(WT_intra_B.shape)
        print(WT_intra_B.head())

        WT_inter_B = WT_loops_anno_addType[WT_loops_anno_addType['type'] == 'inter_B']
        print(WT_inter_B.shape)
        print(WT_inter_B.head())

        WT_inter_mix = WT_loops_anno_addType[WT_loops_anno_addType['type'] == 'inter_mix']
        print(WT_inter_mix.shape)
        print(WT_inter_mix.head())

        WT_other = WT_loops_anno_addType[WT_loops_anno_addType['type'] == 'other']
        print(WT_other.shape)
        print(WT_other.head())

        HS_loops_anno_addType_raw = pd.read_csv(
            src_dir / f'ANNO_CDS0{HS_num}D_combined_clusters_cis_filterChrom_addLoopSpan_filterPETcount{PETcount_thresh}_filterLoopspanNone_addType.txt',
            sep='\t')
        HS_loops_anno_addType = HS_loops_anno_addType_raw[HS_loops_anno_addType_raw['count'] >= scatterPlot_PETthresh]
        print(HS_loops_anno_addType.shape)
        HS_intra_A = HS_loops_anno_addType[HS_loops_anno_addType['type'] == 'intra_A']
        print(HS_intra_A.shape)
        print(HS_intra_A.head())

        HS_inter_A = HS_loops_anno_addType[HS_loops_anno_addType['type'] == 'inter_A']
        print(HS_inter_A.shape)
        print(HS_inter_A.head())

        HS_intra_B = HS_loops_anno_addType[HS_loops_anno_addType['type'] == 'intra_B']
        print(HS_intra_B.shape)
        print(HS_intra_B.head())

        HS_inter_B = HS_loops_anno_addType[HS_loops_anno_addType['type'] == 'inter_B']
        print(HS_inter_B.shape)
        print(HS_inter_B.head())

        HS_inter_mix = HS_loops_anno_addType[HS_loops_anno_addType['type'] == 'inter_mix']
        print(HS_inter_mix.shape)
        print(HS_inter_mix.head())

        HS_other = HS_loops_anno_addType[HS_loops_anno_addType['type'] == 'other']
        print(HS_other.shape)
        print(HS_other.head())
        
        fig, ax = plt.subplots(figsize=(5, 5), dpi=200)
        for spine in ax.spines.values():
            spine.set_linewidth(1)
            spine.set_color('black')

        # Hide the right and top spines
        for spine_position in ['right', 'top']:
            ax.spines[spine_position].set_visible(False)

        # Set the facecolor of the plot
        ax.set_facecolor('white')

        WT_HS_inter_A = [np.log(WT_inter_A['count']), np.log10(HS_inter_A['count'])]
        ax.boxplot(WT_HS_inter_A)

        plt.show()



















