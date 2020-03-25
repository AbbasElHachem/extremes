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

from _00_additional_functions import (select_df_within_period,
                                      resampleDf)
from _08_combine_all_netatmo_into_one_df import max_ppt_thr

register_matplotlib_converters()


plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 12})

plt.ioff()

myFmt = mdates.DateFormatter('%Y-%m-%d')
#==============================================================================
# Read data
#==============================================================================
path_to_daily_netatmo_ppt = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\all_netatmo_ppt_data_daily_.csv"
path_to_daily_dwd_ppt = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\all_dwd_daily_ppt_data_combined_2014_2019_.csv"


# path_to_hourly_netatmo_ppt = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing.csv"
# path_to_hourly_dwd_ppt = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\all_dwd_hourly_ppt_data_combined_2014_2019_.csv"

# path_to_hourly_netatmo_ppt = r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.csv"
# path_to_hourly_dwd_ppt = r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\ppt_dwd_2014_2019_60min_no_freezing_5deg.csv"

path_to_hourly_netatmo_ppt = (
    r'/run/media/abbas/EL Hachem 2019/home_office'
    r'/2020_10_03_Rheinland_Pfalz/ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.csv')
    
path_to_hourly_dwd_ppt = (
    r'/run/media/abbas/EL Hachem 2019/home_office'
    r'/2020_10_03_Rheinland_Pfalz/ppt_dwd_2014_2019_60min_no_freezing_5deg.csv'
    )

# path to stns in RH
path_to_hourly_dwd_coords_in_rh = (
    r'/run/media/abbas/EL Hachem 2019/home_office'
    r'/2020_10_03_Rheinland_Pfalz/dwd_coords_only_in_RH.csv'
    )

# if True, stations how measure high value consecutively
# are investigated and removed if often high values seen
remove_bad_hours_ = True

# select data only within this period (same as netatmo / dwd)
start_date = '2014-01-01 00:00:00'
end_date = '2019-12-31 00:00:00'

resample_data = True
resample_frequencies = ['60min', '120min', '180min',
                        '360min', '720min', '1440min']
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

# dwd in RH coords
df_dwd_coords_in_rh = pd.read_csv(
    path_to_hourly_dwd_coords_in_rh, 
    sep=';', index_col=0)

# if resample_data:
#     for agg_freq in resample_frequencies:
#         print('DWD stations resampling for', agg_freq)
#         df_dwd_resampled = resampleDf(df=df_dwd_hourly,
#                                       agg=agg_freq)
#         # bad_hours = stns_netatmo[stns_netatmo.duplicated()]
#         # good_hours = stns_netatmo[~stns_netatmo.duplicated()]
#
#         dwd_maximum_hrs_dates = df_dwd_resampled.max(
#             axis=1).sort_values()[::-1]
#         dwd_max_100_event = dwd_maximum_hrs_dates[:200].sort_index()
#         dwd_max_100_event = pd.DataFrame(data=np.round(dwd_max_100_event.values, 2),
#                                          index=dwd_max_100_event.index)
#         stns_dwd = df_dwd_resampled.loc[dwd_max_100_event.index, :].idxmax(
#             axis=1)
#         dwd_max_100_event['Station Id'] = stns_dwd.values
#
#     #     netatmo_max_100_event.to_csv(
#     #         r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\netatmo_%s_maximum_100_event.csv" % (
#     #             agg_freq),
#     #         sep=';', header=False)
#     #
#     #     # netatmo_max_100_hours.to_csv(
#     #     #     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\netatmo_hourly_maximum_100_hours_stns.csv",
#     #     #     sep=';', header=False)
#     #
#         dwd_max_100_event.to_csv(
#             r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\dwd_%s_maximum_100_event.csv" % (
#                 agg_freq),
#             sep=';', header=False)
#
#         pass
#
#     raise Exception

# special events for RH

# events = ['29-07-2014', '20-09-2014', '12-06-2015', '17-09-2015',
#           '30-05-2016', '04-06-2016', '03-06-2017', '07-05-2018', '27-05-2018',
#           '11-06-2018', '10-03-2019', '17-03-2019', '21-05-2019', '12-07-2019',
#           '26-07-2019', '27-08-2019']
#
# events_idx = pd.DatetimeIndex(events)
# hourly_range = [pd.date_range(
#     start=(pd.DatetimeIndex([day_start]) - pd.Timedelta(days=1))[0],
#     end=(pd.DatetimeIndex([day_start]) + pd.Timedelta(days=1))[0],
#     freq='60min') for day_start in events_idx]
# hourly_ranges = pd.DataFrame(pd.DatetimeIndex([idx for item in hourly_range
#                                                for idx in item]))
#
# df_dwd_events = df_dwd_hourly.loc[hourly_ranges.values.ravel(), :]
# df_dwd_events2 = df_dwd_events[
#     df_dwd_events.values > 5].drop_duplicates().max(axis=1)
#
# df_dwd_events2 = pd.DataFrame(index=df_dwd_events2.index,
#                               data=df_dwd_events2.values)
# stns_dwd_events = df_dwd_events.loc[df_dwd_events2.index, :].idxmax(
#     axis=1)
#
# # df_dwd_events2['Station ID'] = np.empty(shape=stns_dwd_events.shape[0])
# # for ix in stns_dwd_events.index:
# df_dwd_events2.loc[:, 'Station ID'] = stns_dwd_events.values
# #
# df_dwd_events2.to_csv(
#     r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\dwd_hourly_special_events_5mm_.csv",
#     sep=';', header=False)

#==============================================================================
#
#==============================================================================
# df_netatmo_hourly200 = df_netatmo_hourly[df_netatmo_hourly < 200]
# df_netatmo_hourly200 = df_netatmo_hourly.copy()
# for istn, stn in enumerate(df_netatmo_hourly.columns):
#     for ix, vals in enumerate(df_netatmo_hourly[stn].values):
#         if vals > 200:
#             df_netatmo_hourly200.iloc[ix, istn] = np.nan
# 
# df_netatmo_hourly200.to_csv(path_to_hourly_netatmo_ppt, sep=';', float_format='%0.2f')

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

# select stns in RH
df_dwd_hourly = df_dwd_hourly.loc[
    :,df_dwd_coords_in_rh.index]

dwd_maximum_hrs_dates = df_dwd_hourly.max(axis=1).sort_values()[::-1]
dwd_max_100_hours = dwd_maximum_hrs_dates[:100].sort_index()
dwd_max_100_hours = pd.DataFrame(data=dwd_max_100_hours.values,
                                 index=dwd_max_100_hours.index)
stns_dwd = df_dwd_hourly.loc[dwd_max_100_hours.index, :].idxmax(
    axis=1)
dwd_max_100_hours['Station Id'] = stns_dwd.values


# netatmo_max_100_hours.to_csv(
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\netatmo_%s_maximum_100_event.csv" % (
#         agg_freq),
#     sep=';', header=False)

# netatmo_max_100_hours.to_csv(
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\netatmo_hourly_maximum_100_hours_stns.csv",
#     sep=';', header=False)

# for bw
# dwd_max_100_hours.to_csv(
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\dwd_hourly_maximum_100_hours.csv",
#     sep=';', header=False)
# for RH
# dwd_max_100_hours.to_csv(
#     r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\dwd_hourly_maximum_100_hours.csv",
#     sep=';', header=False)

# dwd_max_100_hours.to_csv(
<<<<<<< HEAD
#     r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\dwd_hourly_maximum_100_hours.csv",
#     sep=';', header=False)

=======
#    r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\dwd_hourly_maximum_100_hours_in_rh.csv",
#    sep=';', header=False)
>>>>>>> branch 'master' of https://github.com/AbbasElHachem/extremes.git
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

# netatmo_max_100_days.to_csv(
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\netatmo_daily_maximum_100_days.csv",
#     sep=';', header=False)
# dwd_max_100_days.to_csv(
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\dwd_daily_maximum_100_days.csv",
#     sep=';', header=False)


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
