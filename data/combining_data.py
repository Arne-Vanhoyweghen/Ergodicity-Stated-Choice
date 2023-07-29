import numpy as np
import pandas as pd
import scipy.io

# Read data from CSV and create a DataFrame with appropriate column names
columns = ['gamble1_1', 'gamble1_2', 'gamble2_1', 'gamble2_2']
gamble_data = pd.read_csv('betcouples.csv', header=None, usecols=[1, 3, 7, 9], names=columns)
mul_data = gamble_data.iloc[40:, :]
linerise_mul_data = mul_data.multiply(1000)
gamble_data.iloc[40:, :] = linerise_mul_data

# Create the 'dynamic' column with 'additive' for the first 40 rows and 'multiplicative' for the rest
gamble_data['wealth'] = 1000
gamble_data['dynamic'] = 0
gamble_data.loc[40:, 'dynamic'] = 1

groups = {'non_timed': {'n_participants': 41, 'nTrials': 80, 'file_name': 'Matrix_Not_Timed'},
          'timed':     {'n_participants': 40, 'nTrials': 80, 'file_name': 'Matrix_Timed'}}

for g in groups:
    df_res = pd.read_excel(f'{groups[g]["file_name"]}.xlsx', header=None)
    df_melted = pd.melt(df_res.T, var_name='participant', value_name='response')
    duplicated_gamble_data = pd.concat([gamble_data] * groups[g]['n_participants'], ignore_index=True)
    merged_df = pd.concat([df_melted, duplicated_gamble_data], axis=1)
    datadict = dict()
    for i in range(groups[g]['n_participants']):
        for c in [0,1]:
            df = merged_df.query('dynamic == @c and participant == @i')
            txt_append = '_add' if c == 0 else '_mul'

            datadict.setdefault(f'x1_1{txt_append}', []).append(np.array(df["gamble1_1"]))
            datadict.setdefault(f'x1_2{txt_append}', []).append(np.array(df["gamble1_2"]))
            datadict.setdefault(f'x2_1{txt_append}', []).append(np.array(df["gamble2_1"]))
            datadict.setdefault(f'x2_2{txt_append}', []).append(np.array(df["gamble2_2"]))
            datadict.setdefault(f'wealth{txt_append}', []).append(np.array(df["wealth"]))
            datadict.setdefault(f'choice{txt_append}', []).append(np.array(df["response"]))
    scipy.io.savemat((f"{g}_data.mat"), datadict, oned_as="row")
