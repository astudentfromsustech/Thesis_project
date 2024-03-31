import pathlib as p
import pandas as pd
import fanc

import numpy as np
# import seaborn as sns
# import matplotlib
# matplotlib.use('Qt5Agg')
# import matplotlib.pyplot as plt
# plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

# def calulate_total_pairs(hic, resolution=1000, chro_list=[]):

def calculate_totalPair_chromList(hic, chro_list):

    def calculate_totalPair_betweenChroms(hic, chro_source, chro_sink):
        totalPair_chroms = 0
        totalEdge_chroms = 0
        for edge in hic.edges((chro_source, chro_sink), lazy=True, norm=False):
            totalPair_chroms += edge.weight
            totalEdge_chroms += 1
        return totalPair_chroms

    totalPair_intra = 0
    totalPair_inter = 0
    n = 0
    for chro_source in chro_list:
        # print(n)
        # print(chro_source)
        for chro_sink in chro_list[n:]:
            # print(chro_source, chro_sink)
            if chro_source == chro_sink:
                totalPair_intra += calculate_totalPair_betweenChroms(hic, chro_source, chro_sink)
                # print(totalPair_intra)
            else:
                totalPair_inter += calculate_totalPair_betweenChroms(hic, chro_source, chro_sink)
                # print(totalPair_inter)
        n += 1
    totalPair = totalPair_intra + totalPair_inter

    return totalPair, totalPair_intra, totalPair_inter

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir_PNAS = root.parent / 'data' / 'PNAS_S2_hic'
    src_dir_home = root.parent / 'data' / 'combined_hic'

    # binSize = 1000
    binSize = 10000


    # hic = fanc.load(src_dir_PNAS / f'WT_combined_PNAS_resolution10000_28000000.hic')
    hic = fanc.load(src_dir_PNAS / f'WT_combined_PNAS.hic')
    print(type(hic))
    # totalPair, totalPair_intra, totalPair_inter = calculate_totalPair_chromList(hic, hic.chromosomes())
    # print(totalPair)
    # print(totalPair_intra)
    # print(totalPair_inter)

    # hic = fanc.load(src_dir_home / f'CDS005D_combined.hic@binSize')
    # hic_dnsample = hic.downsample(0.8)
    # print(type(hic_dnsample))
    # print(hic_dnsample.chromosomes())







