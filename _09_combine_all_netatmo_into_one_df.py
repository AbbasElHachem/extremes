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

from _00_additional_functions import list_all_full_path
from _05_NetatmoData_get_simultaneous_events import (
    split_one_df_file_to_get_stn_id, split_df_file_to_get_alls_stn_ids)

import pandas as pd

dfs_loc = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_bw_grosser_1hour'

dfs_list = list_all_full_path('.csv', dfs_loc)
dfs_list_ppt = list(filter(lambda x: 'coords' not in x, dfs_list))

stn_ids = split_df_file_to_get_alls_stn_ids(dfs_list_ppt)

date_range = pd.date_range('2014-05-01 00:00:00',
                           '2019-05-06 00:00:00',
                           freq='H')
df_all = pd.DataFrame(index=date_range, columns=stn_ids)

for df_file in dfs_list_ppt:
    stn_id = split_one_df_file_to_get_stn_id(df_file)
    print(stn_id)
    in_df = pd.read_csv(df_file,
                        index_col=0,
                        sep=';',
                        parse_dates=True,
                        infer_datetime_format=True,
                        engine='c')
    # FIX ME TO SLOWWWWWW
    df_all.loc[in_df.index, stn_id] = in_df.values.ravel()
#     break
