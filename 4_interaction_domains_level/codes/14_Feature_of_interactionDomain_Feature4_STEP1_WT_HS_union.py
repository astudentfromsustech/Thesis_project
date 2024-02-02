import pathlib as p
import pandas as pd

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
plt.style.use('ggplot')


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def union_intervals(df1, df2):
    # Concatenate the two dataframes
    df = pd.concat([df1[['start', 'end']], df2[['start', 'end']]])

    # Sort by the 'start' column
    df = df.sort_values(by='start')

    # Iterate and merge intervals
    merged = []
    for _, row in df.iterrows():
        if not merged or merged[-1][1] < row['start']:
            # If no overlap, add the interval as is
            merged.append([row['start'], row['end']])
        else:
            # If there is an overlap, merge with the last interval
            merged[-1][1] = max(merged[-1][1], row['end'])

    # Convert to DataFrame
    return pd.DataFrame(merged, columns=['start', 'end'])


if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / '11_call_interaction_domains_STEP3_dropIntervalLen'
    dest_dir = root.parent / 'results' / '14_Feature_of_interactionDomain_Feature4_STEP1_WT_HS_union'
    dest_dir.mkdir(parents=True, exist_ok=True)

    num_list = [['01', '02'], ['04', '05'], ['07', '08']]
    # num_list = [['01', '02']]
    for num_pair in num_list:
        print(num_pair)
        WT_num = num_pair[0]
        HS_num = num_pair[1]


        WT_domain_raw = pd.read_csv(src_dir / f'CDS0{WT_num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
        WT_domain_raw.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        # print(WT_domain_raw.shape)
        # print(WT_domain_raw.head())
        WT_domain = WT_domain_raw[['chromosome', 'start', 'end']]
        # print(WT_domain.shape)
        # print(WT_domain.head())


        HS_domain_raw = pd.read_csv(src_dir / f'CDS0{HS_num}D_filter15P_intervalLen15000.bed', header=None, sep='\t')
        HS_domain_raw.columns = ['chromosome', 'start', 'end', 'accumulated_count']
        # print(HS_domain_raw.shape)
        # print(HS_domain_raw.head())
        HS_domain = HS_domain_raw[['chromosome', 'start', 'end']]
        # print(HS_domain.shape)
        # print(HS_domain.head())

        union = pd.DataFrame()
        chro_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']
        # chro_list = ['chr2L']
        for chro in chro_list:
            print(chro)
            WT_domain_chro = WT_domain[WT_domain['chromosome'] == chro]
            # print(WT_domain_chro.shape)
            # print(WT_domain_chro.head())

            HS_domain_chro = HS_domain[HS_domain['chromosome'] == chro]
            # print(HS_domain_chro.shape)
            # print(HS_domain_chro.head())

            union_chro = union_intervals(WT_domain_chro, HS_domain_chro)
            print(union_chro.shape)
            union_chro.insert(0, 'chromosome', chro)
            union = pd.concat([union, union_chro])
        print(union.shape)
        union['ID'] = [f'{i:03d}' for i in range(1, len(union) + 1)]
        print(union.shape)
        print(union)
        union.to_csv(dest_dir / f'CDS0{WT_num}D_CDS0{HS_num}D_unionIntervals.txt', index=None, sep='\t')


