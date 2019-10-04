#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 08:08:15 2019

@author: abbas
"""

#%%
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pandas.plotting import register_matplotlib_converters

from _00_additional_functions import select_df_within_period
register_matplotlib_converters()


plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 12})

plt.ioff()

myFmt = mdates.DateFormatter('%Y-%m-%d')
#==============================================================================
# Read data
#==============================================================================
path_to_daily_netatmo_ppt = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\all_netatmo_ppt_data_daily_.csv"
path_to_daily_dwd_ppt = r"F:\download_DWD_data_recent\all_dwd_daily_ppt_data_combined_2014_2019_.csv"


path_to_hourly_netatmo_ppt = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing.csv"
path_to_hourly_dwd_ppt = r"F:\download_DWD_data_recent\all_dwd_hourly_ppt_data_combined_2014_2019_.csv"

remove_bad_hours_ = True

# select data only within this period (same as netatmo)
start_date = '2015-01-01 00:00:00'
end_date = '2019-09-30 00:00:00'

#==============================================================================
# # hourly data
#==============================================================================
df_netatmo_hourly = pd.read_csv(
    path_to_hourly_netatmo_ppt,
    sep=';', index_col=0, parse_dates=True,
    infer_datetime_format=True,
    engine='c').dropna(how='all')


df_dwd_hourly = pd.read_csv(
    path_to_hourly_dwd_ppt,
    sep=';', index_col=0, parse_dates=True,
    infer_datetime_format=True,
    engine='c').dropna(how='all')

df_netatmo_hourly = select_df_within_period(df_netatmo_hourly,
                                            start_date,
                                            end_date)

df_dwd_hourly = select_df_within_period(df_dwd_hourly,
                                        start_date,
                                        end_date)
#==============================================================================
#
#==============================================================================
netatmo_maximum_hrs_dates = df_netatmo_hourly.max(axis=1).sort_values()[::-1]
netatmo_max_100_hours = netatmo_maximum_hrs_dates[:100].sort_index()
netatmo_max_100_hours = pd.DataFrame(data=netatmo_max_100_hours.values,
                                     index=netatmo_max_100_hours.index)
stns_netatmo = df_netatmo_hourly.loc[netatmo_max_100_hours.index, :].idxmax(
    axis=1)

netatmo_max_100_hours['Station Id'] = stns_netatmo.values

df_count_event_per_stn = (netatmo_max_100_hours.fillna('')
                          .groupby(netatmo_max_100_hours.columns.tolist()).apply(len)
                          .rename('group_count')
                          .reset_index()
                          .replace('', np.nan)
                          .sort_values(by=['group_count'], ascending=False))

df_many_duplicates_per_stn = df_count_event_per_stn[
    df_count_event_per_stn['group_count'] > 3]
bad_stns = df_many_duplicates_per_stn['Station Id'].values

# bad_hours = stns_netatmo[stns_netatmo.duplicated()]
# good_hours = stns_netatmo[~stns_netatmo.duplicated()]

dwd_maximum_hrs_dates = df_dwd_hourly.max(axis=1).sort_values()[::-1]
dwd_max_100_hours = dwd_maximum_hrs_dates[:100].sort_index()
dwd_max_100_hours = pd.DataFrame(data=dwd_max_100_hours.values,
                                 index=dwd_max_100_hours.index)
stns_dwd = df_dwd_hourly.loc[dwd_max_100_hours.index, :].idxmax(
    axis=1)
dwd_max_100_hours['Station Id'] = stns_dwd.values


netatmo_max_100_hours.to_csv(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\netatmo_hourly_maximum_100_hours.csv",
    sep=';', header=False)

# netatmo_max_100_hours.to_csv(
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\netatmo_hourly_maximum_100_hours_stns.csv",
#     sep=';', header=False)

dwd_max_100_hours.to_csv(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\dwd_hourly_maximum_100_hours.csv",
    sep=';', header=False)

#==============================================================================
# # daily data
#==============================================================================
df_netatmo_daily = pd.read_csv(
    path_to_daily_netatmo_ppt,
    sep=';', index_col=0, parse_dates=True,
    infer_datetime_format=True,
    engine='c').dropna(how='all')
#%%

df_dwd_daily = pd.read_csv(
    path_to_daily_dwd_ppt,
    sep=';', index_col=0, parse_dates=True,
    infer_datetime_format=True,
    engine='c').dropna(how='all')

netatmo_maximum_dates = df_netatmo_daily.max(axis=1).sort_values()[::-1]
dwd_maximum_dates = df_dwd_daily.max(axis=1).sort_values()[::-1]

# dwd_corr_netatmo_max = df_dwd_daily.loc[
#     netatmo_maximum_dates.index, :].max(axis=1).sort_values()[::-1].dropna(how='all')
# netatmo_corr_dwd_max = df_netatmo_daily.loc[
#     dwd_maximum_dates.index, :].max(axis=1).sort_values()[::-1].dropna(how='all')
#
netatmo_max_100_days = netatmo_maximum_dates[:100].sort_index()
# dwd_corr_netatmo_max_100_days = dwd_corr_netatmo_max[:100].sort_index()
#

netatmo_max_100_days = pd.DataFrame(data=netatmo_max_100_days.values,
                                    index=netatmo_max_100_days.index)
stns_netatmo_daily = df_netatmo_daily.loc[netatmo_max_100_days.index, :].idxmax(
    axis=1)
netatmo_max_100_days['Station Id'] = stns_netatmo_daily.values


df_bad_stn_days = netatmo_max_100_days[netatmo_max_100_days[0] > 280]

bad_stns_day_max = df_bad_stn_days[df_bad_stn_days['Station Id'].duplicated()]
df_bad_stn_days = df_bad_stn_days.drop_duplicates(subset=['Station Id'])
df_bad_stn_days = df_bad_stn_days[df_bad_stn_days[0] > 320]

bad_stns_day = df_bad_stn_days['Station Id'].values


dwd_max_100_days = dwd_maximum_dates[:100].sort_index()
# net_corr_dwd_max_100_days = netatmo_corr_dwd_max[:100].sort_index()
dwd_max_100_days = pd.DataFrame(data=dwd_max_100_days.values,
                                index=dwd_max_100_days.index)
stns_dwd_daily = df_dwd_daily.loc[dwd_max_100_days.index, :].idxmax(
    axis=1)
dwd_max_100_days['Station Id'] = stns_dwd_daily.values

netatmo_max_100_days.to_csv(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\netatmo_daily_maximum_100_days.csv",
    sep=';', header=False)
dwd_max_100_days.to_csv(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\dwd_daily_maximum_100_days.csv",
    sep=';', header=False)


#==============================================================================
# FILTERING
#==============================================================================

if remove_bad_hours_:
    print('removing stns hourly', bad_stns)
    for idx, stn_id in zip(netatmo_max_100_hours.index,
                           netatmo_max_100_hours['Station Id'].values):

        if stn_id in bad_stns:
            print(idx, stn_id)
            df_netatmo_hourly.loc[idx, stn_id] = np.nan

    for idxd, stn_idd in zip(netatmo_max_100_days.index,
                             netatmo_max_100_days['Station Id'].values):

        if stn_idd in bad_stns_day:
            print(idxd, stn_idd)
            end = idxd + pd.Timedelta(hours=24)
            time_arr_idx = pd.date_range(start=idxd,
                                         end=end, freq='H')
            bad_hours_in_df = df_netatmo_hourly.index.intersection(
                time_arr_idx)
            df_netatmo_hourly.loc[bad_hours_in_df, stn_idd] = np.nan

    print('Saving Dataframe')

#     print(df_netatmo_hourly.isna().sum())
    df_netatmo_hourly.dropna(how='all', inplace=True)
    df_netatmo_hourly.to_csv(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW',
                                          r'ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.csv'),
                             sep=';', float_format='%.2f')  # temperature humidity ppt

    df_netatmo_hourly.reset_index(inplace=True)
    df_netatmo_hourly.rename(columns={'index': 'Time'}, inplace=True)
    df_netatmo_hourly.to_feather(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW',
                                              r'ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.fk'))
