import pathlib as p

import pandas as pd
import numpy as np
import hicstraw

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'data' / 'combined_hic'

    hic = hicstraw.HiCFile(str(src_dir / "WT_combined_PNAS.hic"))
    print(hic.getChromosomes())
    print(hic.getGenomeID())
    print(hic.getResolutions())

    # result = hicstraw.straw('observed', 'NONE', str(src_dir / "inter.hic"), 'X', 'X', 'BP', 1000000)
    # for i in range(len(result)):
    #     print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))