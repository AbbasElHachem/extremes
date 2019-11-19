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
import os
import timeit
import time

import pandas as pd
import numpy as np
#==============================================================================
#
#==============================================================================

# path_to_DWD_temp_data = (
#     r'x:\exchange\ElHachem\DWD_Reutlingen_Temperature\dwd_reutlingen_temp_2014_2019.csv')

path_to_DWD_temp_data = (
    r"X:\exchange\ElHachem\DWD_Reutlingen_Temperature\all_dwd_daily_temp_data_combined_2014_2019.csv")

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

# def temporal resolution of Netatmo
temp_res = '60min'

# read DWD Temp data
in_df_dwd_temp_data = pd.read_csv(
    path_to_DWD_temp_data,
    engine='c',
    sep=';', index_col=0)

# original
in_df_dwd_temp_data.index = pd.to_datetime(in_df_dwd_temp_data.index,
                                           format='%Y-%m-%d')
# in_df_dwd_temp_data.index = pd.to_datetime(in_df_dwd_temp_data.index,
#                                            format='%d-%m-%Y')


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

# in_df_distance_netatmo_dwd_bw = in_df_distance_netatmo_dwd

netatmo_stn_df_no_freezing_vals = in_df_netatmo_ppt_data.copy()

print('Done getting all data\n')

for netatmo_stn in in_df_netatmo_ppt_data.columns:

    print('Netatmo station is', netatmo_stn)
    netatmo_stn_df = in_df_netatmo_ppt_data.loc[:, netatmo_stn].dropna()

    distance_to_dwd_stns = in_df_distance_netatmo_dwd_bw.loc[netatmo_stn, :]
    min_distance = distance_to_dwd_stns.sort_values()[0]
    dwd_stn_id = distance_to_dwd_stns.sort_values().index[0]

    print('Distance to DWD Temp station', np.round(min_distance, 2))
    dwd_temp_data = in_df_dwd_temp_data.loc[:, dwd_stn_id].drop_duplicates()
#     dwd_temp_data = in_df_dwd_temp_data.drop_duplicates()

    dwd_temp_data_freezing = dwd_temp_data[
        dwd_temp_data < freezing_temp].dropna()

    # freezing days and one day before
    dwd_freezing_days = dwd_temp_data_freezing.index

    print('DWD days with freezing temperatures in Netatmo Ppt',
          dwd_freezing_days.shape[0])

    # generate time period based on netatmo freq for freezing days
    intersect_ix_temp_res = []
    for freez_ix in dwd_freezing_days:
        time_range_freez = pd.date_range(start=freez_ix - pd.Timedelta(days=1),
                                         end=freez_ix + pd.Timedelta(days=1),
                                         freq=temp_res)
        intersect_ix_temp_res.append(time_range_freez)

    # make it a datetime index for intersection with netatmo station
    ix_unique_freez = pd.DatetimeIndex(
        np.unique([ix for ix_ser in intersect_ix_temp_res
                   for ix in ix_ser]))

    # find these days in netatmo ppt data
    # replace these days with nans
    intersect_ix = ix_unique_freez.intersection(netatmo_stn_df.index)

    print('%d / %d events were flagged to nans in Netatmo station'
          % (intersect_ix.size, netatmo_stn_df.size))

    netatmo_stn_df_no_freezing_vals.loc[
        intersect_ix, netatmo_stn_df.name] = np.nan

#==============================================================================
#
#==============================================================================

print('\nSaving new df')

netatmo_stn_df_no_freezing_vals.to_csv(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW'
    r'\ppt_all_netatmo_hourly_stns_combined_new_no_freezing.csv',
    sep=';', float_format='%0.2f')

netatmo_stn_df_no_freezing_vals.reset_index(inplace=True)
netatmo_stn_df_no_freezing_vals.rename(columns={'index': 'Time'}, inplace=True)
netatmo_stn_df_no_freezing_vals.to_feather(
    os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW',
                 r'ppt_all_netatmo_hourly_stns_combined_new_no_freezing.fk'))
STOP = timeit.default_timer()  # Ending time
print(('\n****Done with everything on %s.\nTotal run time was'
       ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
