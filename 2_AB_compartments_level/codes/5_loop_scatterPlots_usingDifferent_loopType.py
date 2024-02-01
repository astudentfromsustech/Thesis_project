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
    dest_dir = root.parent / 'results' / '5_loop_scatterPlots_usingDifferent_loopType'
    dest_dir.mkdir(parents=True, exist_ok=True)


    bin_size = 10000

    PETcount_thresh = 1
    # loop_span_thresh = 500000
    scatterPlot_PETthresh = 2 



    num_list = [['01', '02'], ['04', '05'], ['07', '08']]
    # num_list = [['04', '05']]
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
        
        fig, axes = plt.subplots(2, 4, figsize=(12, 5), dpi=200)

        xticks = [4, 5, 6]
        xticklabels = ['10K', '100K', '1M']
        yticks = [1, 2, 3]
        yticklabels = ['10', '100', '1K']

        for ax in axes.flatten():
            ax.set_facecolor('white')
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)

            for spine in ax.spines.values():
                spine.set_color('black')
                spine.set_linewidth(1)


        axes[0, 0].scatter(np.log10(WT_intra_B['loop_span']), np.log10(WT_intra_B['count']), s=1.5, color='#000033',alpha=0.7, label='WT_intra_B')
        axes[0, 0].scatter(np.log10(WT_inter_B['loop_span']), np.log10(WT_inter_B['count']), s=1.5, color='#0000FF',alpha=0.7, label='WT_inter_B')
        axes[0, 0].legend()
        axes[0, 0].legend(frameon=False)
        # axes[0, 0].set_xticklabels([])

        axes[0,1].scatter(np.log10(WT_intra_A['loop_span']), np.log10(WT_intra_A['count']), s=1.5, color='#808000', alpha=0.7, label='WT_intra_A')
        axes[0, 1].scatter(np.log10(WT_inter_A['loop_span']), np.log10(WT_inter_A['count']), s=1.5, color='#008000',alpha=0.7, label='WT_inter_A')
        axes[0, 1].legend()
        axes[0, 1].legend(frameon=False)
        # axes[0,1].set_xticklabels([])

        axes[0, 2].scatter(np.log10(WT_inter_mix['loop_span']), np.log10(WT_inter_mix['count']), s=1.5, color='#800080',alpha=0.7, label='WT_inter_mix')
        axes[0, 2].legend()
        axes[0, 2].legend(frameon=False)

        axes[0, 3].scatter(np.log10(WT_other['loop_span']), np.log10(WT_other['count']), s=1.5, color='#808080',alpha=0.7, label='WT_other')
        axes[0, 3].legend()
        axes[0, 3].legend(frameon=False)
        # axes[0, 2].set_xticklabels([])
        # axes[0,0].xticks([1000, 10000, 100000, 1000000, 10000000], ['1K', '10K', '100K', '1M', '10M'])
        # axes[0,0].yticks([2, 10, 100, 1000], ['2', '10', '100', '1K'])
        # axes[0,0].xlabel('')
        # axes[0,0].ylabel('')


        axes[1, 0].scatter(np.log10(HS_intra_B['loop_span']), np.log10(HS_intra_B['count']), s=1.5, color='#000033',alpha=0.7, label='HS_intra_B')
        axes[1, 0].scatter(np.log10(HS_inter_B['loop_span']), np.log10(HS_inter_B['count']), s=1.5, color='#0000FF',alpha=0.7, label='HS_inter_B')
        axes[1, 0].legend()
        axes[1, 0].legend(frameon=False)
        # axes[1, 0].set_xticklabels([])

        axes[1,1].scatter(np.log10(HS_intra_A['loop_span']), np.log10(HS_intra_A['count']), s=1.5, color='#808000', alpha=0.7, label='HS_intra_A')
        axes[1, 1].scatter(np.log10(HS_inter_A['loop_span']), np.log10(HS_inter_A['count']), s=1.5, color='#008000',alpha=0.7, label='HS_inter_A')
        axes[1, 1].legend()
        axes[1, 1].legend(frameon=False)
        # axes[1,1].set_xticklabels([])

        axes[1, 2].scatter(np.log10(HS_inter_mix['loop_span']), np.log10(HS_inter_mix['count']), s=1.5, color='#800080',alpha=0.7, label='HS_inter_mix')
        axes[1, 2].legend()
        axes[1, 2].legend(frameon=False)

        axes[1, 3].scatter(np.log10(HS_other['loop_span']), np.log10(HS_other['count']), s=1.5, color='#808080', alpha=0.7, label='HS_other')
        axes[1, 3].legend()
        axes[1, 3].legend(frameon=False)

        # axes[1, 2].set_xticklabels([])

        plt.savefig(dest_dir / f'CDS0{WT_num}_CDS0{HS_num}D_drawScatterPlot_WT_andHS_PETthresh{scatterPlot_PETthresh}.png', dpi=300, bbox_inches='tight')
        # plt.show()



















