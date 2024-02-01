import cooler
import pathlib as p
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
# from hmmlearn import hmm

np.set_printoptions(threshold=np.inf, linewidth=np.inf)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)


def create_chromosome_dataframe(chromosome, num_rows, interval):
    chromosome_data = [chromosome] * num_rows
    start_data = [i * interval for i in range(num_rows)]
    end_data = [(i + 1) * interval - 1 for i in range(num_rows)]

    # Create the DataFrame
    df = pd.DataFrame({
        'chromosome': chromosome_data,
        'start': start_data,
        'end': end_data
    })
    return df

def remove_sparse_rows_and_columns(df):
    # Define the threshold for removing
    threshold = 0.99

    # Calculate the number of zeros in each row and column
    row_zeros = (df == 0).sum(axis=1)
    col_zeros = (df == 0).sum(axis=0)

    # Calculate the total number of entries in each row and column
    total_rows = df.shape[1]  # Number of columns
    total_cols = df.shape[0]  # Number of rows

    # Find rows and columns where the number of zeros is greater than 30% of the entries
    rows_to_remove = row_zeros[row_zeros > threshold * total_rows].index
    cols_to_remove = col_zeros[col_zeros > threshold * total_cols].index

    # Drop the rows and columns from the dataframe
    df_cleaned = df.drop(index=rows_to_remove, columns=cols_to_remove)

    return df_cleaned

def merge_df_with_series_preserving_index(df, series, series_name):
    """
    Merge a DataFrame with a Series based on their indices, preserving the original index of the DataFrame.

    Parameters:
    df (pd.DataFrame): The DataFrame to merge.
    series (pd.Series): The Series to merge.

    Returns:
    pd.DataFrame: Merged DataFrame with the Series as a new column, preserving the original DataFrame index.
    """
    # Preserve the original index of the DataFrame
    original_index = df.index

    # Reset the index of the DataFrame and Series for merging
    df_reset = df.reset_index()
    series_reset = series.reset_index()
    series_reset.columns = ['Index', series_name]  # Rename columns for clarity

    # Merge the DataFrame and Series on the index
    merged_df = pd.merge(df_reset, series_reset, left_on='index', right_on='Index', how='left')

    # Set the original index back to the DataFrame
    merged_df.index = original_index

    # Drop the extra columns
    merged_df.drop(['index', 'Index'], axis=1, inplace=True)

    return merged_df

def add_ratio_columns(df):
    df['WT_ratio'] = np.where(
        df['WT_me3_mean'].notna() & df['WT_ac_mean'].notna(),
        np.log2(df['WT_me3_mean'] / df['WT_ac_mean']),
        np.nan
    )

    # Calculate HS_ratio
    df['HS_ratio'] = np.where(
        df['HS_me3_mean'].notna() & df['HS_ac_mean'].notna(),
        np.log2(df['HS_me3_mean'] / df['HS_ac_mean']),
        np.nan
    )
    return df


def add_compartment_columns(df):
    WT_compartment = []
    HS_compartment = []

    for _, row in df.iterrows():
        if pd.isna(row['WT_ratio']) or pd.isna(row['HS_ratio']):
            WT_compartment.append(np.nan)
            HS_compartment.append(np.nan)
        else:
            WT_compartment.append('A' if row['WT_ratio'] < 0 else 'B')
            HS_compartment.append('A' if row['HS_ratio'] < 0 else 'B')
    df['WT_compartment'] = WT_compartment
    df['HS_compartment'] = HS_compartment
    return df

def add_compartmentCorrected_columns(df):
    WT_compartmentCorrected = []
    HS_compartmentCorrected = []
    WT_Chro_P5 = df['WT_ratio'].abs().quantile(0.05)
    print(WT_Chro_P5)
    HS_Chro_P5 = df['HS_ratio'].abs().quantile(0.05)
    print(HS_Chro_P5)
    for _, row in df.iterrows():
        if pd.isna(row['WT_compartment']) or pd.isna(row['HS_compartment']):
            WT_compartmentCorrected.append(np.nan)
            HS_compartmentCorrected.append(np.nan)
        elif (row['WT_compartment'] != row['HS_compartment']) and ((abs(row['WT_ratio']) <= WT_Chro_P5) or (abs(row['HS_ratio']) <= HS_Chro_P5)):
            if (abs(row['WT_ratio']) <= WT_Chro_P5) and (abs(row['HS_ratio']) <= HS_Chro_P5):
                WT_compartmentCorrected.append(np.nan)
                HS_compartmentCorrected.append(np.nan)
            elif (abs(row['WT_ratio']) <= WT_Chro_P5):
                WT_compartmentCorrected.append(row['HS_compartment'])
                HS_compartmentCorrected.append(row['HS_compartment'])
            else:
                WT_compartmentCorrected.append(row['WT_compartment'])
                HS_compartmentCorrected.append(row['WT_compartment'])
        else:
            WT_compartmentCorrected.append(row['WT_compartment'])
            HS_compartmentCorrected.append(row['HS_compartment'])

    df['WT_compartmentCorrected'] = WT_compartmentCorrected
    df['HS_compartmentCorrected'] = HS_compartmentCorrected
    return df

def add_compartment_switch(df):
    compartment_switch = []
    for _, row in df.iterrows():
        if pd.isna(row['WT_compartmentCorrected']) or pd.isna(row['HS_compartmentCorrected']):
            compartment_switch.append(np.nan)
        elif row['WT_compartmentCorrected'] == row['HS_compartmentCorrected']:
            if row['WT_compartmentCorrected'] == 'A':
                compartment_switch.append('AA')
            else:
                compartment_switch.append('BB')
        elif row['WT_compartmentCorrected'] == 'A':
            compartment_switch.append('AB')
        else:
            compartment_switch.append('BA')
    df['compartment_swich'] = compartment_switch
    return df


def merge_intervals(df):
    # First sort the DataFrame to ensure it's ordered by 'chromosome' and 'start'
    df = df.sort_values(['chromosome', 'start']).reset_index(drop=True)

    # Create an empty DataFrame for merged intervals
    merged_intervals = pd.DataFrame(columns=df.columns)

    # Initialize the first interval
    current_chromosome = df.at[0, 'chromosome']
    current_start = df.at[0, 'start']
    current_end = df.at[0, 'end']

    for index, row in df.iterrows():
        # If the interval is on the same chromosome and the gap between intervals is 1 or less, merge them
        if row['chromosome'] == current_chromosome and row['start'] <= current_end + 1:
            current_end = max(current_end, row['end'])
        else:
            # Add the previous interval
            merged_intervals = pd.concat([merged_intervals,
                                          pd.DataFrame([[current_chromosome, current_start, current_end]],
                                                       columns=df.columns)], ignore_index=True)
            # Start a new interval
            current_chromosome = row['chromosome']
            current_start = row['start']
            current_end = row['end']

    # Add the last interval
    merged_intervals = pd.concat(
        [merged_intervals, pd.DataFrame([[current_chromosome, current_start, current_end]], columns=df.columns)],
        ignore_index=True)
    return merged_intervals

# def plot_cumulative_density(series, title='Cumulative Density Plot', xlabel='Values'):
#     plt.figure(figsize=(8, 6))
#     series.hist(cumulative=True, density=True, bins=100, grid=True)
#     plt.title(title)
#     plt.xlabel(xlabel)
#     plt.ylabel('Cumulative Density')
#     plt.grid(True)
#     plt.show()

if __name__ == '__main__':
    root = p.Path(__file__).absolute().parent

    bin_size = 5000
    dir_num = '5k'

    # bin_size = 10000
    # dir_num = '10k'

    src_dir = root.parent / 'data' / f'combined_cool_{dir_num}_targetChroms'
    dest_dir = root.parent / 'results' / '1_ABcompartments_classification'
    dest_dir.mkdir(parents=True, exist_ok=True)

    num_list = [['01', '04', '02', '05']]
    for num_pair in num_list:
        print(num_pair)
        WT_me3_num = num_pair[0]
        WT_ac_num = num_pair[1]
        print(WT_me3_num)
        print(WT_ac_num)

        HS_me3_num = num_pair[2]
        HS_ac_num = num_pair[3]
        chro_list = ['2L', '2R', '3L', '3R', '4', 'X']
        # chro_list = ['4']

        WT_me3 = cooler.Cooler(str(src_dir / f'CDS0{WT_me3_num}D_combined_{bin_size}_targetChroms.cool'))
        WT_ac = cooler.Cooler(str(src_dir / f'CDS0{WT_ac_num}D_combined_{bin_size}_targetChroms.cool'))

        HS_me3 = cooler.Cooler(str(src_dir / f'CDS0{HS_me3_num}D_combined_{bin_size}_targetChroms.cool'))
        HS_ac = cooler.Cooler(str(src_dir / f'CDS0{HS_ac_num}D_combined_{bin_size}_targetChroms.cool'))

        UNION_df = pd.DataFrame()
        for chro in chro_list:
            print(chro)
            chromosome = 'chr' + chro
            WT_me3_Chro_raw = pd.DataFrame(WT_me3.matrix(balance=False).fetch(chro))
            print(WT_me3_Chro_raw.shape)
            WT_ac_Chro_raw = pd.DataFrame(WT_ac.matrix(balance=False).fetch(chro))
            print(WT_ac_Chro_raw.shape)
            HS_me3_Chro_raw = pd.DataFrame(HS_me3.matrix(balance=False).fetch(chro))
            print(HS_me3_Chro_raw.shape)
            HS_ac_Chro_raw = pd.DataFrame(HS_ac.matrix(balance=False).fetch(chro))
            print(HS_ac_Chro_raw.shape)


            UNION_Chro_df = create_chromosome_dataframe(chromosome, WT_me3_Chro_raw.shape[0], bin_size)
            # print(UNION_Chro_df.shape)
            # print(UNION_Chro_df.head())


            WT_me3_Chro = remove_sparse_rows_and_columns(pd.DataFrame(WT_me3_Chro_raw))
            # print(WT_me3_Chro.shape)
            #
            WT_me3_Chro_normalized = (WT_me3_Chro / WT_me3_Chro.values.mean())
            # print(WT_me3_Chro_normalized.shape)
            # print(WT_me3_Chro_normalized.iloc[0:21, 0:21])
            WT_me3_Chro_mean_colMean = WT_me3_Chro_normalized.mean()
            # print(len(WT_me3_Chro_mean_colMean))
            # print(WT_me3_Chro_mean_colMean)
            #
            WT_ac_Chro = remove_sparse_rows_and_columns(pd.DataFrame(WT_ac_Chro_raw))
            # print(WT_ac_Chro.shape)
            WT_ac_Chro_normalized = (WT_ac_Chro / WT_ac_Chro.values.mean())
            # print(WT_ac_Chro_normalized.shape)
            # print(WT_ac_Chro_normalized.iloc[0:21, 0:21])
            WT_ac_Chro_mean_colMean = WT_ac_Chro_normalized.mean()
            # print(len(WT_ac_Chro_mean_colMean))
            # print(WT_ac_Chro_mean_colMean)
            
            HS_me3_Chro = remove_sparse_rows_and_columns(pd.DataFrame(HS_me3_Chro_raw))
            # print(HS_me3_Chro.shape)
            #
            HS_me3_Chro_normalized = (HS_me3_Chro / HS_me3_Chro.values.mean())
            # print(HS_me3_Chro_normalized.shape)
            # print(HS_me3_Chro_normalized.iloc[0:21, 0:21])
            HS_me3_Chro_mean_colMean = HS_me3_Chro_normalized.mean()
            # print(len(HS_me3_Chro_mean_colMean))
            # print(HS_me3_Chro_mean_colMean)
            #
            HS_ac_Chro = remove_sparse_rows_and_columns(pd.DataFrame(HS_ac_Chro_raw))
            # print(HS_ac_Chro.shape)
            HS_ac_Chro_normalized = (HS_ac_Chro / HS_ac_Chro.values.mean())
            # print(HS_ac_Chro_normalized.shape)
            # print(HS_ac_Chro_normalized.iloc[0:21, 0:21])
            HS_ac_Chro_mean_colMean = HS_ac_Chro_normalized.mean()
            # print(len(HS_ac_Chro_mean_colMean))
            # print(HS_ac_Chro_mean_colMean)
            #

            UNION_Chro_df_addMean = merge_df_with_series_preserving_index(UNION_Chro_df, WT_me3_Chro_mean_colMean, 'WT_me3_mean')
            UNION_Chro_df_addMean = merge_df_with_series_preserving_index(UNION_Chro_df_addMean, WT_ac_Chro_mean_colMean, 'WT_ac_mean')
            UNION_Chro_df_addMean = merge_df_with_series_preserving_index(UNION_Chro_df_addMean, HS_me3_Chro_mean_colMean, 'HS_me3_mean')
            UNION_Chro_df_addMean = merge_df_with_series_preserving_index(UNION_Chro_df_addMean, HS_ac_Chro_mean_colMean, 'HS_ac_mean')
            # print(UNION_Chro_df_addMean.shape)
            # print(UNION_Chro_df_addMean)
            UNION_Chro_df_addMean_addRatio = add_ratio_columns(UNION_Chro_df_addMean)
            # print(UNION_Chro_df_addMean_addRatio.shape)
            # print(UNION_Chro_df_addMean_addRatio)
            UNION_Chro_df_addMean_addRatio_addCompartment = add_compartment_columns(UNION_Chro_df_addMean_addRatio)
            # print(UNION_Chro_df_addMean_addRatio_addCompartment.shape)
            # print(UNION_Chro_df_addMean_addRatio_addCompartment)

            UNION_Chro_df_addMean_addRatio_addCompartment_addCompartmentCorrected = add_compartmentCorrected_columns(UNION_Chro_df_addMean_addRatio_addCompartment)
            # print(UNION_Chro_df_addMean_addRatio_addCompartment_addCompartmentCorrected.shape)
            # print(UNION_Chro_df_addMean_addRatio_addCompartment_addCompartmentCorrected)
            UNION_Chro_df_addMean_addRatio_addCompartment_addCompartmentCorrected_addSwitch = add_compartment_switch(UNION_Chro_df_addMean_addRatio_addCompartment_addCompartmentCorrected)
            # print(UNION_Chro_df_addMean_addRatio_addCompartment_addCompartmentCorrected_addSwitch.shape)
            # print(UNION_Chro_df_addMean_addRatio_addCompartment_addCompartmentCorrected_addSwitch)

            UNION_df = pd.concat([UNION_df, UNION_Chro_df_addMean_addRatio_addCompartment_addCompartmentCorrected_addSwitch])
        print(UNION_df.shape)
        print(UNION_df.head())
        UNION_df.to_csv(dest_dir / f'Union_df_AB_compartment_switch_resolution{bin_size}.txt', index=None, header=True, sep='\t')
        UNION_df_counts = UNION_df['compartment_swich'].value_counts()
        print(UNION_df_counts)

                        # total = sum(UNION_df_counts)
                        # percentages = {k: (v / total) * 100 for k, v in UNION_df_counts.items()}
                        # # Categories as they should appear on the plot from bottom to top
                        # categories = ['AA', 'AB', 'BA', 'BB']
                        # # Colors for each category
                        # colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
                        #
                        # # Creating the plot
                        # fig, ax = plt.subplots()
                        #
                        # # The start location for the percentage chunks
                        # start = 0
                        # for category, color in zip(categories, colors):
                        #     ax.barh('Compartment switches', percentages[category], left=start, color=color, label=category)
                        #     start += percentages[category]
                        # # Adding the legend
                        # ax.legend(loc='lower right', title='Switch Type')
                        # # Setting labels and title
                        # ax.set_xlabel('Percentage of genome under compartment switches')
                        # ax.set_title('Genome Compartment Switch Distribution')
                        # # Display the plot
                        # plt.show()
        print(UNION_df.head(7000))

        WT_score = UNION_df[(UNION_df['WT_compartment'] == 'A') | (UNION_df['WT_compartment'] == 'B')][['chromosome', 'start', 'end', 'WT_ratio', 'WT_compartmentCorrected']]
        HS_score = UNION_df[(UNION_df['HS_compartment'] == 'A') | (UNION_df['HS_compartment'] == 'B')][['chromosome', 'start', 'end', 'HS_ratio', 'HS_compartmentCorrected']]
        print(WT_score.head(30))
        print(HS_score.head(30))
        WT_score.to_csv(dest_dir / f'WT_compartment_score_resolution{bin_size}.bedGraph', header=None, index=None, sep='\t')
        HS_score.to_csv(dest_dir / f'HS_compartment_score_resolution{bin_size}.bedGraph', header=None, index=None, sep='\t')

        ##### produce AB compartment intervals
        WT_A = UNION_df[UNION_df['WT_compartmentCorrected'] == 'A'][['chromosome', 'start', 'end']]
        print(WT_A.shape)
        print(WT_A.head(30))
        WT_A_merged = merge_intervals(WT_A)
        # print(WT_A_merged)
        WT_A_merged.to_csv(dest_dir / f'WT_A_merged_resolution{bin_size}.bed', index=None, header=None, sep='\t')

        HS_A = UNION_df[UNION_df['HS_compartmentCorrected'] == 'A'][['chromosome', 'start', 'end']]
        # print(HS_A.shape)
        # print(HS_A.head(30))
        HS_A_merged = merge_intervals(HS_A)
        # print(HS_A_merged)
        HS_A_merged.to_csv(dest_dir / f'HS_A_merged_resolution{bin_size}.bed', index=None, header=None, sep='\t')

        WT_B = UNION_df[UNION_df['WT_compartmentCorrected'] == 'B'][['chromosome', 'start', 'end']]
        WT_B_merged = merge_intervals(WT_B)
        WT_B_merged.to_csv(dest_dir / f'WT_B_merged_resolution{bin_size}.bed', index=None, header=None, sep='\t')
        HS_B = UNION_df[UNION_df['HS_compartmentCorrected'] == 'B'][['chromosome', 'start', 'end']]
        HS_B_merged = merge_intervals(HS_B)
        HS_B_merged.to_csv(dest_dir / f'HS_B_merged_resolution{bin_size}.bed', index=None, header=None, sep='\t')



