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

netatmo_maximum_hrs_dates = df_netatmo_hourly.max(axis=1).sort_values()[::-1]
netatmo_max_100_hours = netatmo_maximum_hrs_dates[:200].sort_index()
netatmo_max_100_hours = pd.DataFrame(data=netatmo_max_100_hours.values,
                                     index=netatmo_max_100_hours.index)
stns_netatmo = df_netatmo_hourly.loc[netatmo_max_100_hours.index, :].idxmax(
    axis=1)
bad_hours = stns_netatmo[stns_netatmo.duplicated()]
good_hours = stns_netatmo[~stns_netatmo.duplicated()]

netatmo_max_100_hours['Station Id'] = stns_netatmo.values

for ix, stn in enumerate(netatmo_max_100_hours['Station Id'].values):
    if stn in bad_hours.values:
        netatmo_max_100_hours.iloc[ix, 1] = np.nan

netatmo_max_100_hours.dropna(inplace=True)


if remove_bad_hours_:
    #     print(df_netatmo_hourly.isna().sum())
    for idx, stn_id in zip(bad_hours.index, bad_hours.values):
        df_netatmo_hourly.loc[idx, stn_id] = np.nan
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


dwd_maximum_hrs_dates = df_dwd_hourly.max(axis=1).sort_values()[::-1]
dwd_max_100_hours = dwd_maximum_hrs_dates[:200].sort_index()
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
netatmo_max_100_days = netatmo_maximum_dates[:200].sort_index()
# dwd_corr_netatmo_max_100_days = dwd_corr_netatmo_max[:100].sort_index()
#

netatmo_max_100_days = pd.DataFrame(data=netatmo_max_100_days.values,
                                    index=netatmo_max_100_days.index)
stns_netatmo_daily = df_netatmo_daily.loc[netatmo_max_100_days.index, :].idxmax(
    axis=1)
netatmo_max_100_days['Station Id'] = stns_netatmo_daily.values
bad_days = stns_netatmo_daily[stns_netatmo_daily.duplicated()]
good_days = stns_netatmo_daily[~stns_netatmo_daily.duplicated()]

netatmo_max_100_days['Station Id'] = stns_netatmo_daily.values

for ix, stn in enumerate(netatmo_max_100_days['Station Id'].values):
    if stn in bad_days.values:
        netatmo_max_100_days.iloc[ix, 1] = np.nan

netatmo_max_100_days.dropna(inplace=True)


dwd_max_100_days = dwd_maximum_dates[:200].sort_index()
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
# PLOTTING
#==============================================================================

# xticks = pd.date_range(start=netatmo_max_100_days.index[0],
#                        end=netatmo_max_100_days.index[-1], freq='1D')
# fig, ax = plt.subplots(figsize=(30, 16), dpi=150)
# ax.plot(netatmo_max_100_days.index, netatmo_max_100_days.values,
#         label='Netatmo', color='r', alpha=0.75)
#
# ax.plot(dwd_corr_netatmo_max_100_days.index, dwd_corr_netatmo_max_100_days.values,
#         label='DWD', color='b', alpha=0.75)
#
#
# # ax.set_xticks(xticks)
#
#
# ax.set_xlim([netatmo_max_100_days.index[0], netatmo_max_100_days.index[-1]])
#
#
# ax.tick_params(axis='x', rotation=45)
# ax.xaxis.set_major_formatter(myFmt)
# ax.xaxis.set_ticks(
#     xticks[::int(np.round(xticks.shape[0] / 50))])
#
# plt.title('Highest 100 daily maximum rainfall values')
# plt.ylabel('Rainfall mm/d')
# plt.legend(loc=0)
# plt.grid(alpha=0.5)
# plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\daily_maximums.png",
#             frameon=True, papertype='a4',
#             bbox_inches='tight', pad_inches=.2)
