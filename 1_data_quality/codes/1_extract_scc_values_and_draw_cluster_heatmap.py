import pathlib as p
import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use('Qt5Agg')
# matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / "results_hicRep_calculate_scc_pairwise"
    print(src_dir)
    dest_dir = root.parent / 'results' / "1_results_hicRep_scc_cluster_heatmap"
    dest_dir.mkdir(parents=True, exist_ok=True)
    number_samples = 12

    df = pd.DataFrame(index=range(number_samples),columns=range(number_samples))
    df.columns = ['H3K27me3_WT_rep1', 'H3K27me3_HS_rep1', 'H3K27ac_WT_rep1', 'H3K27ac_HS_rep1', 'RNAP2_WT_rep1', 'RNAP2_HS_rep1',
'H3K27me3_WT_rep2', 'H3K27me3_HS_rep2', 'H3K27ac_WT_rep2', 'H3K27ac_HS_rep2', 'RNAP2_WT_rep2', 'RNAP2_HS_rep2']
    df.index = ['H3K27me3_WT_rep1', 'H3K27me3_HS_rep1', 'H3K27ac_WT_rep1', 'H3K27ac_HS_rep1', 'RNAP2_WT_rep1', 'RNAP2_HS_rep1',
'H3K27me3_WT_rep2', 'H3K27me3_HS_rep2', 'H3K27ac_WT_rep2', 'H3K27ac_HS_rep2', 'RNAP2_WT_rep2', 'RNAP2_HS_rep2']
    print(df)
    n = 0
    for file in src_dir.rglob('*.txt'):
        i = n // number_samples
        j = n % number_samples
        print("(i, j): ", i, j)
    #     # print(file)
    #     # scc = pd.read_csv(file, header=None, skiprows=2).iloc[0, 0]
        scc = pd.read_csv(file, header=None, skiprows=2).mean()
        df.iloc[i, j] = scc
        n = n + 1

    print(df)
    # # df = pd.to_numeric(df, errors='coerce')
    # #
    df[['H3K27me3_WT_rep1', 'H3K27me3_HS_rep1', 'H3K27ac_WT_rep1', 'H3K27ac_HS_rep1', 'RNAP2_WT_rep1', 'RNAP2_HS_rep1',
'H3K27me3_WT_rep2', 'H3K27me3_HS_rep2', 'H3K27ac_WT_rep2', 'H3K27ac_HS_rep2', 'RNAP2_WT_rep2', 'RNAP2_HS_rep2']] = df[['H3K27me3_WT_rep1', 'H3K27me3_HS_rep1', 'H3K27ac_WT_rep1', 'H3K27ac_HS_rep1', 'RNAP2_WT_rep1', 'RNAP2_HS_rep1',
'H3K27me3_WT_rep2', 'H3K27me3_HS_rep2', 'H3K27ac_WT_rep2', 'H3K27ac_HS_rep2', 'RNAP2_WT_rep2', 'RNAP2_HS_rep2']].astype('float')
    print(df)

    new_order = ['H3K27me3_WT_rep1', 'H3K27me3_WT_rep2', 'H3K27me3_HS_rep1', 'H3K27me3_HS_rep2', 'H3K27ac_WT_rep1',
                 'H3K27ac_WT_rep2',
                 'H3K27ac_HS_rep1', 'H3K27ac_HS_rep2', 'RNAP2_WT_rep1', 'RNAP2_WT_rep2', 'RNAP2_HS_rep1',
                 'RNAP2_HS_rep2']
    df_reordered = df.loc[new_order, new_order]
    print(df_reordered)

    # ax = sns.clustermap(data=df_reordered,
    #                annot=True,
    #                vmax= 1,
    #                vmin=0,
    #                cmap="vlag",
    #                col_cluster=False,
    #                row_cluster=False
    #                )

    ax = sns.clustermap(data=df_reordered,
                        annot=True,
                        vmax=1,
                        vmin=0,
                        cmap="vlag",
                        col_cluster=False,
                        row_cluster=False,
                        )

    # plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=45)
    ax.ax_heatmap.set_xticks([])
    ax.ax_heatmap.set_yticks([])
    ax.ax_heatmap.set_xticklabels([])
    ax.ax_heatmap.set_yticklabels([])

    ax.cax.set_visible(False)

    # Adjust the heatmap's position to cover the colorbar's space
    heatmap_pos = ax.ax_heatmap.get_position()
    ax.ax_heatmap.set_position([heatmap_pos.x0, heatmap_pos.y0, heatmap_pos.width * 1.05, heatmap_pos.height])
    # plt.show()
    ax.savefig(dest_dir / 'hicRep_scc_cluster_heatmap_noColorBar_noLabel.png', dpi=300)



