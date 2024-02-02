 # import core packages
import warnings
warnings.filterwarnings("ignore")
from itertools import combinations

import pathlib as p
import matplotlib.pyplot as plt
from matplotlib import colors

plt.style.use('seaborn-poster')
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

# import open2c libraries
import bioframe
import cooler
import cooltools
from packaging import version

if version.parse(cooltools.__version__) < version.parse('0.5.2'):
    raise AssertionError("tutorial relies on cooltools version 0.5.2 or higher,"+
    "please check your cooltools version and update to the latest")

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'home_data' / 'combined_cool_10k_targetChroms_normalize0_1_correctionKR_forNormalize0_1'
    annotation_dir = root.parent / 'data' / 'annotations'
    view_df = pd.read_csv(annotation_dir / 'chrom_size_dm3.txt', index_col=0, sep='\t')
    print(view_df)
    resolution = 1000
    num = '01'
    clr = cooler.Cooler(str(src_dir / f'CDS0{num}D_combined_10000_targetChroms_normalize0_1_correctionKR.cool'))
    print(clr.info)

    print(clr.bins()[:10])  # Display the first 10 bins

    # print(clr.pixels()[:10])  # Display the first 10 pixels

    print(clr.chromnames)  # List of chromosome names

    # print(clr.shape)  # Shape of the matrix
    # print(clr_filterChroms.shape)

    # cvd == contacts-vs-distance
    cvd = cooltools.expected_cis(
        clr=clr,
        view_df=view_df,
        smooth=True,
        aggregate_smoothed=True,
        nproc=4  # if you do not have multiple cores available, set to 1
    )
    print(cvd)

    cvd['s_bp'] = cvd['dist'] * resolution
    f, ax = plt.subplots(1, 1)
    for region in view_df['name']:
         ax.loglog(
             cvd['s_bp'].loc[cvd['region1'] == region],
             cvd['balanced.avg'].loc[cvd['region1'] == region],
         )
         ax.set(
             xlabel='separation, bp',
             ylabel='IC contact frequency')
         ax.set_aspect(1.0)
         ax.grid(lw=0.5)
    f.savefig('tets.png', dpi=300, bbox_inches='tight')