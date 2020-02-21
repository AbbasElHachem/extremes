#==============================================================================
# TODO: ADD DOCUMENATION ORGANIZE
#==============================================================================


import datetime
import os
import time
import timeit
import fnmatch
import pandas as pd
import requests
import zipfile
import numpy as np
import math

from _00_additional_functions import (list_all_full_path, resampleDf,
                                      select_df_within_period)

#==============================================================================
#
#==============================================================================
build_one_df = True
delete_df_files = False

make_hdf5_dataset_new = False

temporal_freq = 'hourly'  # '1_minute' '10_minutes' 'hourly' 'daily'
time_period = 'historical'  # 'recent' 'historical'
temp_accr = 'stundenwerte'  # '1minuten' '10minuten' 'stundenwerte'  'tages'

ppt_act = 'RR'  # nieder 'rr' 'RR' 'RR'

stn_name_len = 5  # length of dwd stns Id, ex: 00023

start_date = '2019-01-01 00:00:00'
end_date = '2020-01-30 23:59:00'

temp_freq = 'H'  # '60Min'  # 'Min' '10Min' 'H' 'D'

# '%Y%m%d%H'  # when reading download dwd station df '%Y%m%d%H%M'
time_fmt = '%Y%m%d%H'  # :%S' # for 10Min data use: '%Y%m%d%H%M:%S'

# name of station id and rainfall column name in df
stn_id_col_name = 'STATIONS_ID'
ppt_col_name = 'TT_TU'  # 'TT_TU'  # RS_01 'RWS_10' '  R1' 'RS'
temp_col_name1 = ' LUFTTEMPERATUR'
temp_col_name2 = 'TT_TU'

#==============================================================================
#
#==============================================================================
# main_dir = os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\download_DWD_data_recent')

main_dir = r'F:\dwd'
os.chdir(main_dir)
#==============================================================================
out_dir = r'F:\dwd\pres_hourly_temp'
# if not os.path.exists(out_dir):
#     os.mkdir(out_dir)

os.chdir(out_dir)

# get all df files as list

#     all_df_files = list_all_full_path('.txt', out_extract_df_dir)

all_var_files = list_all_full_path('.txt', out_dir)
all_data_files = [_fi for _fi in all_var_files if 'produkt' in _fi]
all_stn_names = list(set([_fi.split('\\')[-1].split('_')[-1].split('.')[0]
                          for _fi in all_data_files]))
print('total nbr of stations is ', len(all_stn_names))


# create daterange and df full of nans

date_range = pd.date_range(start=start_date, end=end_date, freq=temp_freq)
data = np.zeros(shape=(date_range.shape[0], len(all_stn_names)))
data[data == 0] = np.nan

final_df_combined = pd.DataFrame(data=data,
                                 index=date_range,
                                 columns=all_stn_names)

# read df stn file and fill final df for all stns

for ix, df_txt_file in enumerate(all_data_files):
    print('\n File nbr ', ix, '/', len(all_data_files))

    stn_name = df_txt_file.split('\\')[-1].split('_')[-1].split('.')[0]
    if stn_name in all_stn_names:
        print('\n## Station ID is ##', stn_name)

        # TODO: FIX ME
        in_df = pd.read_csv(df_txt_file, sep=';',
                            dtype=str, engine='c', memory_map=True).dropna()

        try:
            in_df.index = pd.to_datetime(
                in_df[' MESS_DATUM'], format='%Y%m%d%H')
        except KeyError:
            try:
                in_df.index = pd.to_datetime(
                    in_df['MESS_DATUM'], format='%Y%m%d%H')
            except Exception as msg:
                print(msg, 'error with time idx')
                print('check index column')
                raise Exception
        assert int(stn_name) == int(in_df.loc[
            :, stn_id_col_name].values[0].strip(' '))
        assert stn_name in final_df_combined.columns

        in_df = select_df_within_period(in_df, start_date, end_date)

        if in_df.values.size > 0:
            print('\n++\n  Data shape is++',  in_df.values.size)
            try:
                var_data_str = in_df[temp_col_name1].values
                var_data_arr = [np.float(val.strip(' '))
                                for val in var_data_str]
            except Exception as msg:
                pass
                try:
                    var_data_str = in_df[temp_col_name2].values.ravel()
                    var_data_arr = [np.float(val.strip(' '))
                                    for val in var_data_str]
                except Exception as msg2:
                    print('no Var column name found, check again')
                    raise Exception

            try:
                final_df_combined.loc[in_df.index, stn_name] = var_data_arr
            except Exception:
                print('Error')

print('++-- DONE WITH ALL FILES ++--')
# final_df_combined.plot(legend=False)
print('Saving Dataframe')
final_df_combined.dropna(how='all', inplace=True)
# #final_df_combined.set_index('Time', inplace=True, drop=True)
final_df_combined.to_csv(
    os.path.join(main_dir,
                 r'all_dwd_60min_temp_data_combined_2019_2020.csv'),
    sep=';')
# # all_dwd_60min_temp_data_combined_2010_2018.csv
# #     final_df_combined.reset_index(inplace=True)
# #     final_df_combined.rename({'index': 'Time'}, inplace=True)
# # #
# #     final_df_combined.to_feather(
# #         os.path.join(out_dir,
# #                      'all_dwd_10min_ppt_data_combined_2018_2019.fk'))
#
print('done creating df')
# #==============================================================================
# # # In[129]:
# #==============================================================================
#
#
# if make_hdf5_dataset_new:
#     from a_functions import create_hdf5
#
#     print('\n+++\n reading df \n+++\n')
#     final_df_combined = pd.read_csv(
#         os.path.join(
#             main_dir, r'all_dwd_hourly_ppt_data_combined_1995_2019.csv'),
#         #                             r'all_ppt_data_combined.csv'),
#         sep=';',
#         index_col=0)
#     final_df_combined.index = pd.to_datetime(final_df_combined.index,
#                                              format='%Y-%m-%d %H:%M:%S')
#     print('\n+++\n resampling df \n+++\n')
#
#     for temp_freq_resample in freqs_list:
#         if temp_freq_resample == '60Min':
#             df_resampled = final_df_combined
#         else:
#             print(temp_freq_resample)
#             df_resampled = resampleDf(
#                 final_df_combined,
#                 temp_freq_resample)
#     #     ppt_data_all = np.array(df_60min.values)
#         df_resampled = df_resampled[df_resampled >= 0]
#         print('\n+++\n creating HDF5 file \n+++\n')
#
#         str_date = str(df_resampled.index[0]).replace(
#             '-', '').replace(':', '').replace(' ', '')
#         ed_date = str(df_resampled.index[-1]).replace(
#             '-', '').replace(':', '').replace(' ', '')
#
#         output_hdf5_file = r"DWD_%s_ppt_stns_%s_%s_new.h5" % (
#             temp_freq_resample, str_date, ed_date)
#
#         create_hdf5(hf5_file=output_hdf5_file,
#                     start_dt=df_resampled.index[0],
#                     end_dt=df_resampled.index[-1],
#                     nstats=len(df_resampled.columns),
#                     freq=temp_freq_resample,
#                     data_title=(r'Precipitation DWD Data Historical and Recent'
#                                 r'Aggregation %s') % temp_freq_resample,
#                     in_ppt_df=df_resampled,
#                     in_stns_df=stn_names_df,
#                     utm=False)
