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
                                      split_df_file_to_get_alls_stn_ids)

#==============================================================================
#
#==============================================================================
# rain_bw_1hour'
# dfs_loc = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\humidity_bw_1hour'

# dfs_loc = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\rain_bw_1hour'

dfs_loc = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\temperature_bw_1hour'

dfs_list = list_all_full_path('.csv', dfs_loc)
dfs_list_ppt = list(filter(lambda x: 'coords' not in x, dfs_list))

stn_ids = split_df_file_to_get_alls_stn_ids(dfs_list_ppt)
# 2014-04-01 00:00:00 for ppt
date_range = pd.date_range('2012-01-01 00:00:00',
                           '2019-07-31 00:00:00',
                           freq='H')  # 'H'

data_mtx = np.zeros(shape=(date_range.shape[0], len(stn_ids))).astype('float')
data_mtx[data_mtx == 0] = np.nan

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
    in_df.dropna(inplace=True)
    print('Data has the following shape', in_df.values.shape)
    try:
        #         idx_int = in_df.index.intersection(df_all.index)  # .ravel()
        #         df_cmn = in_df.loc[idx_int,].dropna(how='any')
        #         vals_df = in_df.loc[idx_int]  # .values.ravel()
        df_all.loc[in_df.index, stn_id] = in_df.values.ravel()
    except Exception:
        print(stn_id, ' problem in file, skipping file')
        continue

    all_dfs_len -= 1

df_all.to_csv(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW',
                           r'temperature_all_netatmo_hourly_stns_combined_.csv'),
              sep=';')  # temperature humidity ppt

print('done with everything')
