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

import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as md

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from pandas.plotting import register_matplotlib_converters

from _00_additional_functions import (resample_intersect_2_dfs, resampleDf)

from b_get_data import HDF5

rc('font', size=13)
rc('font', family='serif')
rc('axes', labelsize=13)
rcParams['axes.labelpad'] = 13

register_matplotlib_converters()
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
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\disaggregate_DWD_NetAtmo')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)


max_ppt_thr = 100.

#==============================================================================
#
#==============================================================================


def plot_original_disaggregated_values(stn_1_id,
                                       stn_2_id, df_dwd_hourly_orig,
                                       df_dwd_hourly_disagg,
                                       df_netatmo_hourly_orig,
                                       sep_dist_netatmo_dwd,
                                       out_dir):

    fig = plt.figure(figsize=(24, 12), dpi=200)
    ax = fig.add_subplot(111)

    ax.plot(df_dwd_hourly_orig.index,
            df_dwd_hourly_orig.values,
            c='b',
            marker='o',
            # linestyle='--',
            linewidth=2,
            alpha=0.25,
            markersize=3,
            label='DWD Original %s' % stn_1_id)

    ax.plot(df_dwd_hourly_disagg.index,
            df_dwd_hourly_disagg.values,
            c='r',
            marker='+',
            # linestyle='--',
            linewidth=2,
            alpha=0.25,
            markersize=3,
            label='DWD Disaggregated %s' % stn_1_id)

    ax.plot(df_netatmo_hourly_orig.index,
            df_netatmo_hourly_orig.values,
            c='g',
            marker='*',
            # linestyle='--',
            linewidth=2,
            alpha=0.25,
            markersize=3,
            label='Netatmo Original %s' % stn_2_id)

    xfmt = md.DateFormatter('%Y-%m-%d')

    ax.xaxis.set_major_locator(MultipleLocator(25))
    ax.yaxis.set_major_locator(MultipleLocator(2))

    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(0.0, max(df_dwd_hourly_orig.values.max(),
                         df_dwd_hourly_disagg.values.max(),
                         df_netatmo_hourly_orig.values.max()) + 1)

    ax.set_ylabel('Stn  %s   Precipitation  in mm/hour ' % (stn_1_id))
    ax.tick_params('y', colors='darkblue')

    ax.set_title("Stn: %s vs Stn: %s; \n Distance: %0.1f m; "
                 % (stn_1_id, stn_2_id,
                     sep_dist_netatmo_dwd))
    ax.legend(loc=0)
    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.5)

    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,
                             'orig_vs_disaggregated_ppt_stn_%s_vs_stn_%s_.png'
                             % (stn_1_id, stn_2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure ')
    return

#==============================================================================
#
#==============================================================================


def scatter_original_disaggregated_values(stn_1_id,
                                          stn_2_id,
                                          df_dwd_hourly_orig,
                                          df_dwd_hourly_disagg,
                                          df_netatmo_hourly_orig,
                                          sep_dist_netatmo_dwd,
                                          out_dir):
    plt.ioff()
    fig = plt.figure(figsize=(24, 12), dpi=200)
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
#     ax2.plot(time_arr, df2.values, c='red', marker='+', markersize=2,
#             alpha=0.25)  # , label=stn2_id)

    ax.scatter(df_dwd_hourly_disagg.values.ravel(),
               df_dwd_hourly_orig.values.ravel(),
               c='r',
               marker='+',
               # linestyle='--',
               # linewidth=2,
               alpha=0.5,
               s=8,
               label='DWD Original vs Disaggregated %s' % stn_1_id)

    ax2.scatter(df_dwd_hourly_disagg.values.ravel(),
                df_netatmo_hourly_orig.values.ravel(),
                c='g',
                marker='*',
                # linestyle='--',
                # linewidth=2,
                alpha=0.25,
                s=5,
                label='Netatmo Original %s vs DWD Disaggregated %s'
                % (stn_2_id, stn_1_id))

    # calculate correlations (pearson and spearman)
    pear_corr_dwd_orig_dwd_disagg = pears(df_dwd_hourly_disagg.values,
                                          df_dwd_hourly_orig.values)[0]

    spr_corr_dwd_orig_dwd_disagg = spr(df_dwd_hourly_disagg.values,
                                       df_dwd_hourly_orig.values)[0]

    pear_corr_dwd_orig_netatmo = pears(df_netatmo_hourly_orig.values,
                                       df_dwd_hourly_orig.values)[0]

    spr_corr_dwd_orig_netatmo = spr(df_netatmo_hourly_orig.values,
                                    df_dwd_hourly_orig.values)[0]

    pear_corr_netatmo_dwd_disagg = pears(df_dwd_hourly_disagg.values,
                                         df_netatmo_hourly_orig.values)[0]

    spr_corr_netatmo_dwd_disagg = spr(df_dwd_hourly_disagg.values,
                                      df_netatmo_hourly_orig.values)[0]

    # plot 45 deg line
    _min = min(0, min(df_dwd_hourly_orig.values.min(),
                      df_dwd_hourly_disagg.values.min(),
                      df_netatmo_hourly_orig.values.min()))

    _max = max(df_dwd_hourly_orig.values.max(),
               df_dwd_hourly_disagg.values.max(),
               df_netatmo_hourly_orig.values.max())

    ax.plot([_min, _max], [_min, _max], c='k', linestyle='--', alpha=0.25)
    ax2.plot([_min, _max], [_min, _max], c='b', linestyle='--', alpha=0.25)

    # set plot limit
    ax.set_xlim(-0.1, _max + 0.1)
    ax.set_ylim(-0.1, _max + 0.1)
    ax2.set_ylim(-0.1, _max + 0.1)

    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_major_locator(MultipleLocator(2))

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_xlabel('Stn  %s  Disaggregated Precipitation  in mm/hour '
                  % (stn_1_id))
    ax.set_ylabel('Stn  %s  Orig DWD Ppt in mm/hour ' % (stn_1_id))
    ax2.set_ylabel('Stn  %s  Netatmo Ppt in mm/hour ' % (stn_2_id))

    ax.tick_params('y', colors='red')
    ax2.tick_params('y', colors='g')

    ax.set_title("Hourly Stn: %s vs Stn: %s; Distance: %0.1f m \n "
                 " DWD Orig- DWD Disagg:"
                 " Pearson Corr %0.2f ; Spearman Corr %0.2f \n"

                 " DWD Orig- Netatmo:"
                 " Pearson Corr %0.2f ; Spearman Corr %0.2f \n"

                 " Netatmo- DWD Disagg:"
                 " Pearson Corr %0.2f ; Spearman Corr %0.2f"
                 % (stn_1_id, stn_2_id,
                     sep_dist_netatmo_dwd,
                     pear_corr_dwd_orig_dwd_disagg,
                    spr_corr_dwd_orig_dwd_disagg,
                    pear_corr_dwd_orig_netatmo,
                    spr_corr_dwd_orig_netatmo,
                    pear_corr_netatmo_dwd_disagg,
                    spr_corr_netatmo_dwd_disagg))
#     ax.legend(loc=0)
#     ax2.legend(loc=0)
    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.5)

#     plt.xticks(rotation=45)

#     plt.tight_layout()
    plt.savefig(
        os.path.join(out_dir,
                     'scatter_orig_vs_disaggregated_ppt_stn_%s_vs_stn_%s_.png'
                     % (stn_1_id, stn_2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure ')
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

                # df_combined = pd.DataFrame(columns=[stn_id, stn_2_id])

                df_netatmo_hourly, df_dwd_hourly = resample_intersect_2_dfs(
                    idf1,
                    idf2,
                    '60min')
#                 df_dwd_hourly['Time'] = df_dwd_hourly.index
                if len(df_netatmo_hourly.values) > 0:
                    #                     raise Exception

                    df_dwd_daily = resampleDf(df_dwd_hourly, 'D')
                    df_netatmo_daily = resampleDf(df_netatmo_hourly, 'D')

                    empty_data_arr = np.empty(shape=(df_dwd_hourly.shape[0]))
                    empty_data_arr[empty_data_arr == 0] = np.nan

                    df_disaggregated = pd.DataFrame(data=empty_data_arr,
                                                    index=df_dwd_hourly.index)

                    for idx, val in zip(df_dwd_hourly.index,
                                        df_dwd_hourly.values):
                        daily_idx = datetime.date(idx.year, idx.month, idx.day)
                        if val == 0:
                            #                             print('Val is 0')
                            df_disaggregated.loc[idx] = val

                        if val != 0:
                            if daily_idx in df_dwd_daily.index:
                                #                                 print('val is', val, 'at ', daily_idx)
                                sum_daily_dwd = df_dwd_daily.loc[daily_idx].values
                                sum_daily_netatmo = df_netatmo_daily.loc[daily_idx].values
                                # TODO FIX ME
#                                 print('sum dwd is ', sum_daily_dwd)
#                                 print('sum netatmo is', sum_daily_netatmo)
                                if sum_daily_dwd == 0:
                                    ppt_disagg = 0
                                elif sum_daily_netatmo == 0:
                                    ppt_disagg = 0
                                else:
                                    ppt_disagg = (
                                        val / sum_daily_dwd) * sum_daily_netatmo

#                                 print('disaggregated ppt is', ppt_disagg)
                                df_disaggregated.loc[idx] = ppt_disagg
                            else:
                                df_disaggregated.loc[idx] = np.nan
                                print('index not in daily values ppt is nan')

                    plot_original_disaggregated_values(stn_id,
                                                       stn_2_id,
                                                       df_dwd_hourly,
                                                       df_disaggregated,
                                                       df_netatmo_hourly,
                                                       min_dist,
                                                       out_save_dir_orig)
                    scatter_original_disaggregated_values(stn_id,
                                                          stn_2_id,
                                                          df_dwd_hourly,
                                                          df_disaggregated,
                                                          df_netatmo_hourly,
                                                          min_dist,
                                                          out_save_dir_orig)
#                             break
#                     df_combined = pd.DataFrame(
#                         index=df_dwd_hourly.index,
#                         data=df_dwd_hourly.values,
#                         columns=['DWD original'])
#                     df_combined['DWD disaggregated'] = df_disaggregated.values
#                     df_combined['Netatmo original'] = df_netatmo_hourly.values
#                     df_combined.to_csv(
#                         os.path.join(
#                             out_save_dir_orig,
#                             'combined_orig_disaggregated_DWD_%s_Netatmo_%s_hourly_data.csv'
#                             % (stn_id, stn_2_id)),
#                         sep=';', float_format='%.2f')

                else:
                    print('empty df')
                    continue

        except Exception as msg:
            print(msg)


if __name__ == '__main__':
    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    construct_netatmo_dwd_daily_dfs(path_to_ppt_netatmo_data,
                                    path_to_ppt_hdf_data,
                                    distance_matrix_df_file)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
