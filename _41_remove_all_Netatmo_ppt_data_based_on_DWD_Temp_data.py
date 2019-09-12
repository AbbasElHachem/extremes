# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Read DWD Temperature data in BW and remove all Netatmo station
where at nearby DWD station the temperature was below 5degrees

Created on: 2019-09-12

"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# ===========================================================================

import timeit
import time

import pandas as pd
import numpy as np
#==============================================================================
#
#==============================================================================

path_to_DWD_temp_data = (
    r'F:\DWD_download_Temperature_data'
    r'\all_dwd_daily_temp_data_combined_2014_2019.csv')

path_to_Netatmo_ppt_data = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW'
    r'\ppt_all_netatmo_hourly_stns_combined_new.csv')

path_to_distance_mtx_Netatmo_ppt_DWD_temp = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\distance_mtx_in_m_NetAtmo_ppt_DWD_temp.csv")

freezing_temp = 5  # below 5 degrees it is snow

#==============================================================================
#
#==============================================================================


print('**** Started on %s ****\n' % time.asctime())
START = timeit.default_timer()  # to get the runtime of the program

print('Reading data\n')

# read DWD Temp data
in_df_dwd_temp_data = pd.read_csv(
    path_to_DWD_temp_data,
    engine='c',
    sep=';', index_col=0)

in_df_dwd_temp_data.index = pd.to_datetime(in_df_dwd_temp_data.index,
                                           format='%Y-%m-%d')

# read Netatmo ppt data
in_df_netatmo_ppt_data = pd.read_csv(
    path_to_Netatmo_ppt_data,
    sep=';', index_col=0,
    engine='c')
in_df_netatmo_ppt_data.index = pd.to_datetime(in_df_netatmo_ppt_data.index,
                                              format='%Y-%m-%d %H:%M:%S')

# read distance matrix dwd-netamot ppt
in_df_distance_netatmo_dwd = pd.read_csv(
    path_to_distance_mtx_Netatmo_ppt_DWD_temp,
    sep=';', index_col=0)

# keep only stations in BW with data from (2014-2019)
in_df_distance_netatmo_dwd_bw = in_df_distance_netatmo_dwd.loc[
    in_df_netatmo_ppt_data.columns, in_df_dwd_temp_data.columns]

netatmo_stn_df_no_freezing_vals = in_df_netatmo_ppt_data.copy()

print('Done getting all data\n')

for netatmo_stn in in_df_netatmo_ppt_data.columns:
    print('Netatmo station is', netatmo_stn)
    netatmo_stn_df = in_df_netatmo_ppt_data.loc[:, netatmo_stn].dropna()

    distance_to_dwd_stns = in_df_distance_netatmo_dwd_bw.loc[netatmo_stn, :]
    min_distance = distance_to_dwd_stns.sort_values()[0]
    dwd_stn_id = distance_to_dwd_stns.sort_values().index[0]
    print('Distance to DWD Temp station', np.round(min_distance, 2))
    dwd_ppt_data = in_df_dwd_temp_data.loc[:, dwd_stn_id].drop_duplicates()

    dwd_ppt_data_freezing = dwd_ppt_data[dwd_ppt_data < freezing_temp].dropna()

    dwd_freezing_days = dwd_ppt_data_freezing.index
    dwd_freezing_days_dayb4 = dwd_freezing_days - pd.Timedelta(days=1)

    intersect_ix = dwd_freezing_days.intersection(netatmo_stn_df.index)

    print('DWD days with freezing temperatures in Netatmo Ppt',
          intersect_ix.shape[0])

    intersect_ix_dayb4 = dwd_freezing_days_dayb4.intersection(
        netatmo_stn_df.index)

    netatmo_stn_df_no_freezing_vals.loc[
        intersect_ix, netatmo_stn_df.name] = np.nan

    netatmo_stn_df_no_freezing_vals.loc[
        intersect_ix_dayb4, netatmo_stn_df.name] = np.nan


#     for time_ix in netatmo_stn_df.index:
#
#         ix_year, ix_month, ix_day = time_ix.year, time_ix.month, time_ix.day
#
#         time_ix_short_fmt = pd.datetime(
#             year=ix_year, month=ix_month, day=ix_day)
#
#         ix = np.unique(dwd_ppt_data.index[
#             np.logical_and(dwd_ppt_data.index.year == ix_year,
#                            np.logical_and(dwd_ppt_data.index.month == ix_month,
#                                           dwd_ppt_data.index.day == ix_day))])
#         ix_before = ix - pd.Timedelta(days=1)
#
#         # netatmo_stn_df
#         dwd_temp_val_same_day = dwd_ppt_data.loc[ix].values
#         try:
#             dwd_temp_val_day_before = dwd_ppt_data.loc[ix_before].values
#         except KeyError:
#             dwd_temp_val_day_before = np.nan
#
#         if (dwd_temp_val_same_day < freezing_temp) or (
#                 dwd_temp_val_day_before < freezing_temp):
#             print(time_ix, 'DWD Temp', dwd_temp_val_same_day,
#                   'Netatmo PPT', np.round(netatmo_stn_df.loc[time_ix], 2),
#                   'changed to nan')
#             netatmo_stn_df_no_freezing_vals.loc[
#                 time_ix, netatmo_stn_df.name] = np.nan

print('saving new df')

netatmo_stn_df_no_freezing_vals.to_csv(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW'
    r'\ppt_all_netatmo_hourly_stns_combined_new_no_freezing.csv',
    sep=';', float_format='%0.2f')


STOP = timeit.default_timer()  # Ending time
print(('\n****Done with everything on %s.\nTotal run time was'
       ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
