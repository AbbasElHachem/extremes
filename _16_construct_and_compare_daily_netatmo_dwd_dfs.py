# !/usr/bin/env python.
# -*- coding: utf-8 -*-


"""
Look at pairs of stations (NetAtmo-DWD) 

Select nearest stations, for one DWD station, aggregate to daily
and disaggregate according to relative values of NetAtmo station

Check if this gives reasonable Temporal structure or not 

Save the daily DWD and NetAtmo sums in a dataframe, use it to find 
scaling factor between the values
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"
#==============================================================================
#
#==============================================================================
import os
import timeit
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from _00_additional_functions import (resample_intersect_2_dfs, resampleDf)
from b_get_data import HDF5

plt.ioff()
#==============================================================================
#
#==============================================================================

path_to_ppt_netatmo_data = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                            r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_ppt_netatmo_data), 'wrong NETATMO Ppt file'

path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')
assert os.path.exists(path_to_ppt_hdf_data), 'wrong NETATMO Ppt file'

coords_dwd_df_file = (r'X:\exchange\ElHachem\niederschlag_deutschland'
                      r'\1993_2016_5min_merge_nan.csv')
assert os.path.exists(coords_dwd_df_file), 'wrong DWD coords file'

coords_netatmo_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                          r'\NetAtmo_BW\rain_bw_1hour\netatmo_bw_1hour_coords.csv')

assert os.path.exists(coords_netatmo_df_file), 'wrong NETATMO coords file'


distance_matrix_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                           r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD.csv')
assert os.path.exists(distance_matrix_df_file), 'wrong Distance MTX  file'

out_save_dir_orig = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\compare_plot_cycles_DWD_NetAtmo')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)


max_ppt_thr = 100.

save_df_common = False
#==============================================================================
#
#==============================================================================


def plot_values_2_stns(stn1_id, stn2_id, df1, df2, seperate_dist,
                       temp_freq, out_dir):

    assert df1.shape[0] == df2.shape[0], 'Dfs have different shapes'
    plt.ioff()
    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    time_vals = df1.index.to_pydatetime()
    time_arr = md.date2num(time_vals)

    ax.plot(time_arr, df1.values, c='r', marker='+', markersize=2,
            alpha=0.5, label=stn1_id)
    ax.plot(time_arr, df2.values, c='g', marker='*', markersize=2,
            alpha=0.5, label=stn2_id)
    xfmt = md.DateFormatter('%Y-%m-%d')

    ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.yaxis.set_major_locator(MultipleLocator(3))

    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(-0.1, max(df1.values.ravel().max(),
                          df2.values.ravel().max()) + 1)
    ax.set_ylabel('Ppt mm/day ',
                  rotation=-90)
    ax.tick_params('y', colors='darkred')

    ax.set_title("Ratio between values Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; " % (stn1_id, stn2_id,
                                      seperate_dist,
                                      temp_freq))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.xticks(rotation=30)
    ax.legend(loc='upper right')
#     ax2.legend(loc='best')
#     plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_plot_values_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Line plot')

    return

#==============================================================================
#
#==============================================================================


def plot_ratio_2_stns(stn1_id, stn2_id, df1, df2, seperate_dist,
                      temp_freq, out_dir):

    assert df1.shape[0] == df2.shape[0], 'Dfs have different shapes'
    plt.ioff()
    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    time_vals = df1.index.to_pydatetime()
    time_arr = md.date2num(time_vals)
    ratio_values = df1.values.ravel() / df2.values.ravel()
    ratio_values = np.nan_to_num(ratio_values)
    ratio_values[ratio_values > 1000] = 0
#     ratio_values[(0 <= ratio_values) & (ratio_values < 100)]

    ax.plot(time_arr, ratio_values, c='darkblue', marker='o', markersize=2,
            alpha=0.25)  # , label=stn1_id)

    xfmt = md.DateFormatter('%Y-%m-%d')

    ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.yaxis.set_major_locator(MultipleLocator(3))

    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(- 0.01, ratio_values.max() + 1)
#     ax.set_xlabel('Time')
    ax.set_ylabel(r'Ratio daily values Stn  %s  Stn  %s ' %
                  (stn1_id, stn2_id))

    ax.tick_params('y', colors='darkblue')

    ax.set_title("Ratio between values Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; " % (stn1_id, stn2_id,
                                      seperate_dist,
                                      temp_freq))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.xticks(rotation=30)

#     ax2.legend(loc='best')
#     plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_ratio_values_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Line plot')

    return

#==============================================================================
#
#==============================================================================


def plot_max_daily_vals_2_stns(stn1_id, stn2_id, df1, df2, seperate_dist,
                               temp_freq, out_dir):

    assert df1.shape[0] == df2.shape[0], 'Dfs have different shapes'

    ppt_df_netatmo = df1.groupby([df1.index.dayofyear]).max()
    ppt_df_dwd = df2.groupby([df2.index.dayofyear]).max()

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    time_vals = ppt_df_netatmo.index

#     ratio_values[(0 <= ratio_values) & (ratio_values < 100)]

    ax.plot(time_vals, ppt_df_netatmo.values, c='darkblue',
            marker='o', markersize=2,
            alpha=0.5, label=stn1_id)

    ax.plot(time_vals, ppt_df_dwd.values, c='darkred',
            marker='+', markersize=2,
            alpha=0.5, label=stn2_id)
#     xfmt = md.DateFormatter('%Y-%m-%d')

    ax.xaxis.set_major_locator(MultipleLocator(15))
    ax.yaxis.set_major_locator(MultipleLocator(1))

#     ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(-0.1, max(ppt_df_netatmo.values.ravel().max(),
                          ppt_df_dwd.values.ravel().max()) + 1)
    ax.set_xlim(0, 365)
#     ax.set_xlabel('Time')
    ax.set_ylabel(r'Average Ppt (mm/day) ')
    ax.set_xlabel('Day of year')
    ax.tick_params('y', colors='darkblue')

    ax.set_title("Maximum daily value Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; \n Start: %s - End: %s"
                 % (stn1_id, stn2_id,
                     seperate_dist,
                     temp_freq,
                     str(df1.index[0]),
                     str(df1.index[-1])))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.xticks(rotation=30)

    ax.legend(loc='best')
#     plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_max_daily_values_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Line plot')

    return

#==============================================================================
#
#==============================================================================


def plot_max_monthly_vals_2_stns(stn1_id, stn2_id, df1, df2, seperate_dist,
                                 temp_freq, out_dir):

    assert df1.shape[0] == df2.shape[0], 'Dfs have different shapes'
    ppt_netatmo = df1.resample('M',
                               label='right',
                               closed='right').sum()
#     netatmo_monthly_max = ppt_netatmo.groupby([ppt_netatmo.index.month]).max()
    netatmo_monthly_avg = ppt_netatmo.groupby([ppt_netatmo.index.month]).mean()
#     netatmo_monthly_min = ppt_netatmo.groupby([ppt_netatmo.index.month]).min()

    ppt_df_dwd = df2.resample('M',
                              label='right',
                              closed='right').sum()
#     dwd_monthly_max = ppt_df_dwd.groupby([ppt_df_dwd.index.month]).max()
    dwd_monthly_avg = ppt_df_dwd.groupby([ppt_df_dwd.index.month]).mean()
#     dwd_monthly_min = ppt_df_dwd.groupby([ppt_df_dwd.index.month]).min()
    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

#     time_vals = ppt_df_netatmo.index.to_pydatetime()
#     time_arr = md.date2num(time_vals)

#     ax.plot(netatmo_monthly_max.index,
#             netatmo_monthly_max.values, c='darkred',
#             marker='o', markersize=4,
#             alpha=0.5, label=stn1_id)
#
#     ax.plot(dwd_monthly_max.index, dwd_monthly_max.values, c='darkblue',
#             marker='o', markersize=4,
#             alpha=0.5, label=stn2_id)

    ax.plot(netatmo_monthly_avg.index,
            netatmo_monthly_avg.values, c='r',
            marker='+', markersize=4,
            alpha=0.5)  # , label=stn1_id)

    ax.plot(dwd_monthly_avg.index, dwd_monthly_avg.values, c='b',
            marker='+', markersize=4,
            alpha=0.5)  # , label=stn2_id)

#     ax.plot(netatmo_monthly_min.index,
#             netatmo_monthly_min.values, c='salmon',
#             marker='*', markersize=4,
#             alpha=0.5)  # , label=stn1_id)
#
#     ax.plot(dwd_monthly_min.index, dwd_monthly_min.values, c='lightblue',
#             marker='*', markersize=4,
#             alpha=0.5)  # , label=stn2_id)

#     xfmt = md.DateFormatter('%d-%m-%y')
#
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(5))
#
#     ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(-0.1, max(netatmo_monthly_avg.values.ravel().max(),
                          dwd_monthly_avg.values.ravel().max()) + 1)
#     ax.set_ylim(-0.1, max(netatmo_monthly_max.values.ravel().max(),
#                           dwd_monthly_max.values.ravel().max()) + 1)
#     ax.set_xlim(0, 365)
#     ax.set_xlabel('Time')
    ax.set_ylabel(r'Ppt (mm/month) ')
    ax.set_xlabel('Month of year')
    ax.tick_params('y', colors='k')

    ax.set_title("Mean monthly value Stn: %s vs Stn: %s;"
                 "\n Distance: %0.1f m; "
                 "Time Freq: %s; \n Start: %s - End: %s"
                 % (stn1_id, stn2_id,
                     seperate_dist,
                     temp_freq,
                     str(df1.index[0]),
                     str(df1.index[-1])))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.xticks(rotation=0)

    ax.legend(loc='best')
#     plt.tight_layout()
    plt.savefig(os.path.join(out_dir,
                             '%s_mean_monthly_values_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Line plot')

    return
#==============================================================================
#
#==============================================================================


def construct_netatmo_dwd_daily_dfs(netatmo_ppt_df_file,
                                    path_to_ppt_hdf_data,
                                    distance_matrix_df_file):
    HDF52 = HDF5(infile=path_to_ppt_hdf_data)

    in_netatmo_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                     index_col=0, sep=';',
                                     parse_dates=True,
                                     infer_datetime_format=True,
                                     engine='c')
    netatmo_stns_ids = in_netatmo_stns_df.columns

    in_df_distance_netatmo_dwd = pd.read_csv(distance_matrix_df_file,
                                             sep=';', index_col=0)

    for stn_id in netatmo_stns_ids:

        print('First Netatmo Stn Id is', stn_id)

        try:
            idf1 = in_netatmo_stns_df.loc[:, stn_id]
            idf1.dropna(axis=0, inplace=True)
            idf1 = idf1[idf1 < max_ppt_thr]

            distances_to_stn1 = in_df_distance_netatmo_dwd.loc[stn_id, :]
            sorted_distances = distances_to_stn1.sort_values(
                ascending=True)

            min_dist = sorted_distances.values[0]
            if min_dist <= 5000:

                stn_2_id = sorted_distances.index[0]

                idf2 = HDF52.get_pandas_dataframe(ids=[stn_2_id])
                idf2 = idf2[idf2 < max_ppt_thr]

                print('Second DWD Stn Id is', stn_2_id)

                df_combined = pd.DataFrame(columns=[stn_id, stn_2_id])

                df_netatmo_hourly, df_dwd_hourly = resample_intersect_2_dfs(
                    idf1,
                    idf2,
                    '60min')

                if (df_netatmo_hourly.values.shape[0] > 1000):
                    if save_df_common:
                        df_combined.loc[:, stn_id] = df_netatmo_hourly
                        df_combined.loc[:, stn_2_id] = df_dwd_hourly

                        df_combined_daily = resampleDf(df_combined, '1440min')

                        df_combined_daily.to_csv(
                            os.path.join(
                                out_save_dir_orig,
                                'combined_daily_df_netatmo_%s_dwd_%s_sep_distance_%.1fm_.csv'
                                % (stn_id, stn_2_id, min_dist)),
                            sep=';', float_format='%0.2f')

                    df_netatmo_daily = resampleDf(df_netatmo_hourly, '1440min')
                    df_dwd_daily = resampleDf(df_dwd_hourly, '1440min')

                    plot_values_2_stns(stn_id,
                                       stn_2_id,
                                       df_netatmo_daily,
                                       df_dwd_daily,
                                       min_dist,
                                       '1440min',
                                       out_save_dir_orig)

                    plot_ratio_2_stns(stn_id,
                                      stn_2_id,
                                      df_netatmo_daily,
                                      df_dwd_daily,
                                      min_dist,
                                      '1440min',
                                      out_save_dir_orig)

#                     plot_max_daily_vals_2_stns(stn_id,
#                                                stn_2_id,
#                                                df_netatmo_daily,
#                                                df_dwd_daily,
#                                                min_dist,
#                                                '1440min',
#                                                out_save_dir_orig)

                    plot_max_monthly_vals_2_stns(stn_id,
                                                 stn_2_id,
                                                 df_netatmo_daily,
                                                 df_dwd_daily,
                                                 min_dist,
                                                 '1440min',
                                                 out_save_dir_orig)

                else:
                    print('empty df')
                    continue

        except Exception as msg:
            print(msg)
#         break


if __name__ == '__main__':
    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    construct_netatmo_dwd_daily_dfs(path_to_ppt_netatmo_data,
                                    path_to_ppt_hdf_data,
                                    distance_matrix_df_file)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
