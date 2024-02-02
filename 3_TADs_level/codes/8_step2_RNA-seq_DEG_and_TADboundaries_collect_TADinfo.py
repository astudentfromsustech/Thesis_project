import pathlib as p
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 500)

def determine_label(row, intervals):
    tss = row['TSS']
    for _, interval in intervals.iterrows():
        if interval['Start'] <= tss <= interval['End']:
            return interval['Label']

def calculate_promoter(row, promoter_length=1000):
    if row['strand'] == '+':
        return row['start'] - promoter_length, row['start']
    else:
        return row['end'], row['end'] + promoter_length
def skip_first_row(group):
    return group.iloc[1:]

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent
    src_dir = root.parent / 'results' / 'home_data_results'
    dest_dir = root.parent / 'results' / 'home_data_results' / '8_step2_RNA-seq_DEG_and_TADboundaries_collect_TADinfo'
    dest_dir.mkdir(parents=True, exist_ok=True)
    thresh = 0.05
    delta = 0.01
    
    num_list = ['01', '02', '04', '05', '07', '08']
    # num_list = ['01']
    for num in num_list:
        print(num)
    #num = '01'

        TADs_raw = pd.read_csv(src_dir / f'CDS0{num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh{thresh}_delta{delta}_dfr' /
                       f'CDS0{num}D_domains.bed', header=None, sep='\t')
        TADs_raw.columns = ['chromosome', 'start', 'end', 'ID', 'score_left', 'dot', 'start2', 'end2', 'color']
        #print(TADs_raw.shape)
        #print(TADs_raw.head())
        TADs = TADs_raw[['chromosome', 'start', 'end', 'score_left']]
        print(TADs.shape)
        print(TADs.head())

        sep_scores = pd.read_csv(src_dir / f'CDS0{num}D_combined_10000_targetChroms_normalize0_1_correctionKR_min30k_max100k_step10k_thresh{thresh}_delta{delta}_dfr' /
                       f'CDS0{num}D_boundaries.bed', header=None, sep='\t')
        sep_scores.columns = ['chromosome', 'start', 'end', 'ID', 'score', 'dot']
        # print(sep_scores.shape)
        # print(sep_scores.head())
    
        sep_scores_filterFirstRow = sep_scores.groupby('chromosome').apply(skip_first_row)
        sep_scores_filterFirstRow.reset_index(drop=True, inplace=True)
        # print(sep_scores_filterFirstRow.shape)
        # print(sep_scores_filterFirstRow.head())

        TADs_addScoreRight = TADs.copy()
        TADs_addScoreRight['score_right'] = sep_scores_filterFirstRow['score']
        TADs_addScoreRight['score_right'] = TADs_addScoreRight['score_right'].round(2)
        TADs_addScoreRight['score_left'] = TADs_addScoreRight['score_left'].round(2)
        # print(TADs_addScoreRight.shape)
        # print(TADs_addScoreRight.head())
    
        TADs_addScoreRight['ID'] =  TADs.groupby('chromosome').cumcount() + 1
        print(TADs_addScoreRight.shape)
        print(TADs_addScoreRight.head())

        TADs_addScoreRight['bound_left_start'] = TADs['start'] - 5000
        TADs_addScoreRight['bound_left_end'] = TADs['start'] + 5000
        TADs_addScoreRight['bound_right_start'] = TADs['end'] - 5000
        TADs_addScoreRight['bound_right_end'] = TADs['end'] + 5000
        print(TADs_addScoreRight.shape)
        print(TADs_addScoreRight.head(5))

        TADs_addScoreRight = TADs_addScoreRight[['chromosome', 'start', 'end', 'ID', 'bound_left_start', 'bound_left_end', 'score_left',
                                                 'bound_right_start', 'bound_right_end', 'score_right']]
        print(TADs_addScoreRight.shape)
        print(TADs_addScoreRight.head(5))



        replace_values = {'2L': 'chr2L', '2R': 'chr2R', '3L': 'chr3L', '3R': 'chr3R', '4': 'chr4', 'X': 'chrX'}
        TADs_addScoreRight['chromosome'] = TADs_addScoreRight['chromosome'].replace(replace_values)
        print(TADs_addScoreRight.shape)
        print(TADs_addScoreRight.head())

        TADs_addScoreRight.to_csv(dest_dir / f'CDS0{num}D_combined_allTADsInfo_includingBothSepScore_thresh{thresh}_delta{delta}.txt', index=None, sep='\t')



