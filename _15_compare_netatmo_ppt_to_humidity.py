# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Look at pairs of stations (NetAtmo-NetAtmo) 

Compare Precipitation values to Humidity values of nearest station
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

import numpy as np
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import matplotlib.dates as md


from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears


from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from pandas.plotting import register_matplotlib_converters


from _00_additional_functions import (
    resampleDf, resample_Humidity_Df, constrcut_contingency_table)

register_matplotlib_converters()

plt.ioff()

rc('font', size=13)
rc('font', family='serif')
rc('axes', labelsize=13)
rcParams['axes.labelpad'] = 13

majorLocator = MultipleLocator(2)
minorLocator = MultipleLocator(1)
#==============================================================================
#
#==============================================================================
path_to_ppt_netatmo_data = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                            r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_ppt_netatmo_data), 'wrong NETATMO Ppt file'

ppt_coords_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                      r'\NetAtmo_BW\rain_bw_1hour\netatmo_bw_1hour_coords.csv')

assert os.path.exists(ppt_coords_df_file), 'wrong NETATMO coords file'


path_to_hum_netatmo_data = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                            r'\NetAtmo_BW\humidity_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_hum_netatmo_data), 'wrong NETATMO Hum file'


distance_matrix_df_file_ppt_hum = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\distance_mtx_in_m_NetAtmo_ppt_Netatmo_hum.csv"
assert os.path.exists(distance_matrix_df_file_ppt_hum), 'wrong distance df'

out_save_dir_orig = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_humidity')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

x_col_name = ' lon'
y_col_name = ' lat'


# threshold for CDF, consider only above thr, below is P0
ppt_thr = .5
max_ppt_thr = 100.

ppt_thrs_list = [0.5, 1, 2, 0.5, 5]
hum_thrs_list = [80, 85, 90, 95, 95]

# till 1 day '5min', '10min', '15min', '30min',
aggregation_frequencies = ['60min', '90min', '120min', '180min', '240min',
                           '360min', '480min', '720min', '1440min']

#==============================================================================
#
#==============================================================================


def plt_bar_plot_ppt_hum_stns(stn1_id, stn2_id, seperate_distance,
                              df1, df2, temp_freq, out_dir):
    ''' plot line plots between two stations'''
    print('plotting line plots')

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    time_vals = df1.index.to_pydatetime()
    time_arr = md.date2num(time_vals)

    ax.plot(time_arr, df1.values, c='darkblue', marker='o', markersize=2,
            alpha=0.25)  # , label=stn1_id)

    ax2 = ax.twinx()
    ax2.plot(time_arr, df2.values, c='red', marker='+', markersize=2,
             alpha=0.25)  # , label=stn2_id)

    xfmt = md.DateFormatter('%Y-%m-%d')

    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(10))

    ax2.yaxis.set_major_locator(MultipleLocator(10))

    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(0.0, df1.values.max() + 1)
#     ax.set_xlabel('Time')
    ax.set_ylabel('Stn  %s   Precipitation  in mm/%s ' % (stn1_id, temp_freq))
    ax.tick_params('y', colors='darkblue')

    ax2.set_ylim(0.0, 100.01)
    ax2.set_ylabel('Stn %s   Humidity  in percentage %s'
                   % (stn2_id, temp_freq), rotation=-90)
    ax2.tick_params('y', colors='red')

    ax.set_title("Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; " % (stn1_id, stn2_id,
                                      seperate_distance,
                                      temp_freq))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.xticks(rotation=30)
#     ax.legend(loc='best')
#     ax2.legend(loc='best')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_mean_lineplot_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Line plot')

    return


#==============================================================================
#
#==============================================================================


def plt_bar_plot_ppt_mean_min_max_hum_stns(stn1_id, stn2_id,
                                           seperate_distance,
                                           df1, df2_mean, df2_min, df2_max,
                                           temp_freq, out_dir):
    ''' plot line plots between two stations'''
    print('plotting line plots')

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    time_vals = df1.index.to_pydatetime()
    time_arr = md.date2num(time_vals)

    ax.plot(time_arr, df1.values, c='darkblue', marker='o', markersize=2,
            alpha=0.25)  # , label=stn1_id)

    ax2 = ax.twinx()
    ax2.plot(time_arr, df2_mean.values, c='orange', marker='+', markersize=2,
             alpha=0.5)  # , label=stn2_id)
    ax2.plot(time_arr, df2_min.values, c='g', marker='*', markersize=2,
             alpha=0.5)  # , label=stn2_id)
    ax2.plot(time_arr, df2_max.values, c='r', marker='d', markersize=2,
             alpha=0.5)  # , label=stn2_id)

    xfmt = md.DateFormatter('%Y-%m-%d')

    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(10))

    ax2.yaxis.set_major_locator(MultipleLocator(10))

    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(0.0, df1.values.max() + 1)
#     ax.set_xlabel('Time')
    ax.set_ylabel('Stn  %s   Precipitation  in mm/%s ' % (stn1_id, temp_freq))
    ax.tick_params('y', colors='darkblue')
    ax2.set_ylim(0.0, 100.01)
    ax2.set_ylabel('Stn %s   Humidity  in percentage %s'
                   % (stn2_id, temp_freq), rotation=-90)
    ax2.tick_params('y', colors='darkred')

    ax.set_title("Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; " % (stn1_id, stn2_id,
                                      seperate_distance,
                                      temp_freq))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.xticks(rotation=30)
#     ax.legend(loc='best')
#     ax2.legend(loc='best')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,
                             'mean_%s_all_lineplot_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Line plot')

    return


#==============================================================================
#
#==============================================================================
def plt_scatter_plot_ppt_hum(stn1_id, stn2_id, seperate_distance,
                             df1, df2, ppt_min_thr, temp_freq, out_dir):
    ''' plot scatter plots between two stations'''
    print('plotting scatter plots')

    ppt_abv_thr1 = df1[df1 >= ppt_min_thr]
    ppt_abv_thr2 = df2[df2 >= ppt_min_thr]

    ppt_abv_thr1.dropna(inplace=True)
    ppt_abv_thr2.dropna(inplace=True)

    idx_cmn = ppt_abv_thr1.index.intersection(ppt_abv_thr2.index)

    df1_cmn_df2 = ppt_abv_thr1.loc[idx_cmn]
    df2_cmn_df1 = ppt_abv_thr2.loc[idx_cmn]

    values_x = df1_cmn_df2.values.ravel()
    values_y = df2_cmn_df1.values.ravel()

    # calculate correlations (pearson and spearman)
    corr = pears(values_x, values_y)[0]
    rho = spr(values_x, values_y)[0]

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    ax.scatter(values_x, values_y, c='r', marker='.', alpha=0.35)

    # plot 45 deg line
    _min = min(values_x.min(), values_y.min())
    _max = max(values_x.max(), values_y.max())
    ax.plot([_min, _max], [_min, _max], c='k', linestyle='--', alpha=0.4)

    # set plot limit
    ax.set_xlim(-0.2, _max + 0.01)
    ax.set_ylim(-0.2, _max + 0.01)

    ax.xaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax.set_ylabel('Observed %s Humidity in Percentage station Id: %s'
                  % (temp_freq, stn2_id))
    ax.set_xlabel('Observed %s Rainfall in mm station Id: %s'
                  % (temp_freq, stn1_id))
    ax.set_title("Stn: %s vs Stn: %s; \n Distance: %0.1f m; "
                 "Time Freq: %s; \n Pearson Cor=%0.3f; "
                 "Spearman Cor=%0.3f" % (stn1_id, stn2_id,
                                         seperate_distance,
                                         temp_freq, corr, rho))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,
                             'mean_%s_scatter_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Scatter')

    return


#==============================================================================
#
#==============================================================================

def plot_contingency_tables_as_a_sequence_ppt_hum_stns(stn_id,
                                                       stn_2_id,
                                                       min_dist,
                                                       ppt_thrs_list,
                                                       hum_thrs_list,
                                                       df_stn1,
                                                       df_stn2,
                                                       aggregation_frequencies_lst,
                                                       out_dir):
    '''
    Do contingenccy tables a function of time and thresholds for 2 stations
    '''
    print('Plotting contingency_table a sequence')

    plt.ioff()

    #colors = ['r', 'b', 'g', 'orange', 'k']

    if (df_stn1.values.shape[0] > 0 and
            df_stn2.values.shape[0] > 0):

        for _, (ppt_thr, hum_thr) in enumerate(zip(ppt_thrs_list,
                                                   hum_thrs_list)):

            print('Ppt Threshold is', ppt_thr)

            for temp_freq in aggregation_frequencies_lst:
                fig = plt.figure(figsize=(16, 12), dpi=200)
                ax = fig.add_subplot(111)
                ax.set_aspect(1)
                print('Time freq is', temp_freq)
                df_resampled_ppt = resampleDf(df_stn1, temp_freq)

#                 df_resampled2_hum_mean = resample_Humidity_Df(df_stn2,
#                                                           temp_freq,
#                                                           method='mean')

#                 df_resampled2_hum_min = resample_Humidity_Df(df_stn2,
#                                                              temp_freq,
#                                                              method='min')

                df_resampled2_hum_max = resample_Humidity_Df(df_stn2,
                                                             temp_freq,
                                                             method='max')

                idx_cmn = df_resampled_ppt.index.intersection(
                    df_resampled2_hum_max.index)

                df_ppt_common1 = df_resampled_ppt.loc[idx_cmn]

#                 df_common2_mean = df_resampled2_hum_mean.loc[idx_cmn]
#                     df_common2_min = df_resampled2_hum_min.loc[idx_cmn]
                df_hum_common2_max = df_resampled2_hum_max.loc[idx_cmn]

#                 df_ppt_common1.fillna(0, inplace=True)
#                 df_hum_common2_max.fillna(0, inplace=True)

                if (df_ppt_common1.values.shape[0] > 10 and
                        df_hum_common2_max.values.shape[0] > 10):
                    print('getting percentages for contingency table')

                    df_ppt_common1_new = pd.DataFrame(
                        data=df_ppt_common1.values,
                        index=df_ppt_common1.index,
                        columns=[stn_id])

                    df_hum_common2_max_new = pd.DataFrame(
                        data=df_hum_common2_max.values,
                        index=df_hum_common2_max.index,
                        columns=[stn_2_id])

                    try:
                        (df_both_below_thr,
                         df_first_abv_second_below_thr,
                         df_first_below_second_abv_thr,
                         df_both_abv_thr) = constrcut_contingency_table(stn_id,
                                                                        stn_2_id,
                                                                        df_ppt_common1_new,
                                                                        df_hum_common2_max_new,
                                                                        ppt_thr,
                                                                        hum_thr)

                        conf_arr = np.array([[df_both_below_thr,
                                              df_first_abv_second_below_thr],
                                             [df_first_below_second_abv_thr,
                                              df_both_abv_thr]])
                        print(conf_arr)
                    except Exception as msg:
                        print('error while calculating contingency_table',
                              msg, temp_freq)
                        continue
                else:
                    print('empty df values')
                    break

                if conf_arr.shape[0] > 0:
                    print('plotting for Ppt thr, temp frequency is', temp_freq)

                    sn.heatmap(conf_arr, annot=True, vmin=0.0,
                               vmax=100.0, fmt='.2f')

                    ax.set_title('contingency table \n Humidiy thr %0.0f \n Rainfall %0.1f '
                                 '%s vs %s \n distance: %0.1f m, time freq: %s'
                                 % (hum_thr, ppt_thr, stn_id,
                                    stn_2_id, min_dist, temp_freq))
                    plt.savefig(
                        os.path.join(
                            out_dir,
                            'mean_contingency_table_as_a_sequence_%s_%s_%s_ppt_%0.0f_hum_%0.0f.png'
                            % (stn_id, stn_2_id, temp_freq, ppt_thr, hum_thr)),
                        bbox_inches='tight')
                    plt.close(fig)
            print('Done plotting contingency_tableas a sequence')
        else:
            print('empty df not plotting Contingency table')
    return


#==============================================================================
#
#==============================================================================

def compare_ppt_to_humidity_stns(netatmo_ppt_df_file,
                                 netatmo_humidity_df_file,
                                 distance_matrix_df_file):
    in_netatmo_ppt_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                         index_col=0, sep=';',
                                         parse_dates=True,
                                         infer_datetime_format=True,
                                         engine='c')
    stns_ppt_ids = in_netatmo_ppt_stns_df.columns

    in_netatmo_humidity_stns_df = pd.read_csv(netatmo_humidity_df_file,
                                              index_col=0, sep=';',
                                              parse_dates=True,
                                              infer_datetime_format=True,
                                              engine='c')

    in_df_distance_netatmo_netatmo = pd.read_csv(distance_matrix_df_file,
                                                 sep=';', index_col=0)

    for ppt_stn_id in stns_ppt_ids[5:]:
        print('First Ppt Stn Id is', ppt_stn_id)
        try:
            idf1 = in_netatmo_ppt_stns_df.loc[:, ppt_stn_id]
            idf1.dropna(axis=0, inplace=True)
            idf1 = idf1[idf1 < max_ppt_thr]

            distances_to_stn1 = in_df_distance_netatmo_netatmo.loc[ppt_stn_id, :]
            sorted_distances = distances_to_stn1.sort_values(ascending=True)

            min_dist = sorted_distances.values[0]
            if min_dist <= 5000:

                stn_2_id = sorted_distances.index[0]

                idf2 = in_netatmo_humidity_stns_df.loc[:, stn_2_id]
                idf2.dropna(axis=0, inplace=True)
                print('Second Stn Id is', stn_2_id)
                for tem_freq in aggregation_frequencies:
                    print('Aggregation is: ', tem_freq)

                    df_resampled1 = resampleDf(idf1, tem_freq)
                    df_resampled2_mean = resample_Humidity_Df(idf2,
                                                              tem_freq,
                                                              method='mean')

#                     df_resampled2_min = resample_Humidity_Df(idf2,
#                                                              tem_freq,
#                                                              method='min')

#                     df_resampled2_max = resample_Humidity_Df(idf2,
#                                                              tem_freq,
#                                                              method='max')

                    idx_cmn = df_resampled1.index.intersection(
                        df_resampled2_mean.index)
#
                    df_common1 = df_resampled1.loc[idx_cmn]
                    df_common2_mean = df_resampled2_mean.loc[idx_cmn]

#                     df_common2_min = df_resampled2_min.loc[idx_cmn]
#                     df_common2_max = df_resampled2_max.loc[idx_cmn]
#                     df_common1, df_common2 = resample_intersect_2_dfs(
#                         idf1, idf2, tem_freq)
                    if (df_common1.values.shape[0] > 0 and
                            df_common2_mean.values.shape[0] > 0):

                        out_save_dir = os.path.join(
                            out_save_dir_orig,
                            'ppt_%s_hum_%s' % (ppt_stn_id, stn_2_id))

                        if not os.path.exists(out_save_dir):
                            os.mkdir(out_save_dir)
                        try:
                            pass
                            plt_bar_plot_ppt_hum_stns(ppt_stn_id,
                                                      stn_2_id,
                                                      min_dist,
                                                      df_common1,
                                                      df_common2_mean,
                                                      tem_freq,
                                                      out_save_dir)

#                             plt_bar_plot_ppt_mean_min_max_hum_stns(
#                                 ppt_stn_id,
#                                 stn_2_id,
#                                 min_dist,
#                                 df_common1,
#                                 df_common2_mean,
#                                 df_common2_min,
#                                 df_common2_max,
#                                 tem_freq,
#                                 out_save_dir)

                            plt_scatter_plot_ppt_hum(ppt_stn_id,
                                                     stn_2_id,
                                                     min_dist,
                                                     df_common1,
                                                     df_common2_mean,
                                                     0,
                                                     tem_freq,
                                                     out_save_dir)

                        except Exception as msg:
                            print('error while plotting', msg, tem_freq)
                            continue
                        # break
                    else:
                        print('empty df')
                        continue
                if os.path.exists(out_save_dir):
                    print('out directory is ', out_save_dir)
                    plot_contingency_tables_as_a_sequence_ppt_hum_stns(ppt_stn_id,
                                                                       stn_2_id,
                                                                       min_dist,
                                                                       ppt_thrs_list,
                                                                       hum_thrs_list,
                                                                       idf1,
                                                                       idf2,
                                                                       aggregation_frequencies,
                                                                       out_save_dir)
                else:
                    continue
        except Exception as msg:
            print(msg)

#         break


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    compare_ppt_to_humidity_stns(path_to_ppt_netatmo_data,
                                 path_to_hum_netatmo_data,
                                 distance_matrix_df_file_ppt_hum)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
