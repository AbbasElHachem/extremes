# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Combine all Netatmo data into one Dataframe

Time, Ids,
     ppt
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

#==============================================================================
#
#==============================================================================
import pandas as pd
import os
import numpy as np


from _00_additional_functions import (list_all_full_path,
                                      split_one_df_file_to_get_stn_id,
                                      split_df_file_to_get_alls_stn_ids,
                                      select_df_within_period)

#==============================================================================
#
#==============================================================================
# rain_bw_1hour'
# dfs_loc = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\humidity_bw_1hour'

# dfs_loc = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\rain_bw_1hour'
dfs_loc = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\rain_UK_1hour'
# dfs_loc = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\temperature_bw_1hour'

dfs_list = list_all_full_path('.csv', dfs_loc)
dfs_list_ppt = list(filter(lambda x: 'coords' not in x, dfs_list))

stn_ids = split_df_file_to_get_alls_stn_ids(dfs_list_ppt)
# 2014-04-01 00:00:00 for ppt

list_years = ['2015', '2016', '2017', '2018', '2019']

for _year in list_years:
    date_range = pd.date_range('%s-01-01 00:00:00' % _year,
                               '%s-12-30 00:00:00' % _year,
                               freq='H')  # 'H'

    max_ppt_thr = 100  # maximum ppt values per hour
    initial_vals_to_remove = 4  # most likely calibration or test values
    minimal_number_of_vals = 30 * 24  # 1 month of hourly data

    data_mtx = np.zeros(
        shape=(date_range.shape[0], len(stn_ids))).astype('float')
    data_mtx[data_mtx == 0] = -9  # np.nan

    df_all = pd.DataFrame(index=date_range, columns=stn_ids, data=data_mtx)

    all_dfs_len = len(dfs_list_ppt)

    for df_file in dfs_list_ppt:
        print('Number of files is', all_dfs_len)
        stn_id = split_one_df_file_to_get_stn_id(df_file)
        print(stn_id)
        in_df = pd.read_csv(df_file,
                            index_col=0,
                            sep=';',
                            parse_dates=True,
                            infer_datetime_format=True,
                            engine='c')
        # in_df.dropna(inplace=True)
        in_df = in_df[(0 <= in_df) & (in_df <= max_ppt_thr)]  # (0 <= in_df) &

        in_df = in_df.iloc[initial_vals_to_remove:]  # remove first values
        in_df = select_df_within_period(df=in_df, start=date_range[0],
                                        end=date_range[-1])
        in_df.dropna(inplace=True)
        if in_df.values.shape[0] >= minimal_number_of_vals:
            print('Data has the following shape', in_df.values.shape)
            try:
                #         idx_int = in_df.index.intersection(df_all.index)  # .ravel()
                #         df_cmn = in_df.loc[idx_int,].dropna(how='any')
                #         vals_df = in_df.loc[idx_int]  # .values.ravel()
                df_all.loc[in_df.index, stn_id] = in_df.values.ravel()
            except Exception:
                print(stn_id, ' problem in file, skipping file')
                continue
        else:
            print('Station did not have enough data')

        all_dfs_len -= 1

    print('Saving Dataframe')
    df_all.dropna(how='all', inplace=True)
    df_all.to_csv(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW',
                               r'ppt_all_netatmo_uk_hourly_stns_combined_%s.csv'
                               % _year),  # new
                  sep=';')  # , float_format='%.2f')  # temperature humidity ppt

#     df_all.reset_index(inplace=True)
#     df_all.rename(columns={'index': 'Time'}, inplace=True)
#     df_all.to_feather(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW',
#                                    r'ppt_all_netatmo_hourly_stns_combined_new.fk'))

    print('done with everything')
