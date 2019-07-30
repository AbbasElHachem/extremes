# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Look at pairs of stations (DWD-DWD) 

For different aggregations (15min, 30min, 60min, ..., 6hours, ..., 1 day)
Calculate P0 (probability that rainfall < 1mm)
Construct Cdf for extremes ( > 1mm) and compare
Calculate the correlation between the ranks on all scales
Find if the similarities or differences in the stations are due to 
a systematic or a random error.
If a systematic error is present (bias) find the correction factor 
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
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as md
import seaborn as sn
import scipy.spatial as spatial

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from pandas.plotting import register_matplotlib_converters

from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from collections import OrderedDict

from _00_additional_functions import (resample_intersect_2_dfs,
                                      get_cdf_part_abv_thr,
                                      calculate_probab_ppt_below_thr,
                                      constrcut_contingency_table)

from b_get_data import HDF5

register_matplotlib_converters()

#==============================================================================
#
#==============================================================================
rc('font', size=13)
rc('font', family='serif')
rc('axes', labelsize=13)
rcParams['axes.labelpad'] = 13

majorLocator = MultipleLocator(2)
minorLocator = MultipleLocator(1)

path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')

coords_df_file = (r'X:\exchange\ElHachem\niederschlag_deutschland'
                  r'\1993_2016_5min_merge_nan.csv')
assert os.path.exists(coords_df_file), 'wrong DWD coords file'

out_save_dir_orig = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_DWD_DWD')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

x_col_name = 'Rechtswert'
y_col_name = 'Hochwert'

# threshold for CDF, consider only above thr, below is P0
ppt_thr_min = .5

# used for P0 calculation
ppt_thrs_list = [0.5, 1, 2]


# till 1 day
aggregation_frequencies = ['5min', '10min', '15min', '30min', '60min', '90min',
                           '120min', '180min', '240min',  '360min',
                           '480min', '720min', '1440min']

# which neighbor to chose ,starts with 1
neighbor_to_chose = 1

#==============================================================================
#
#==============================================================================


def get_dwd_stns_coords(coords_df_file, x_col_name, y_col_name, index_col,
                        sep_type):
    '''function used to return to coordinates dataframe of the DWD stations'''
    in_coords_df = pd.read_csv(coords_df_file, sep=sep_type,
                               index_col=index_col, engine='c')
    stn_ids = in_coords_df.index
    x_vals = in_coords_df[x_col_name].values.ravel()
    y_vals = in_coords_df[y_col_name].values.ravel()
    return in_coords_df, x_vals, y_vals, stn_ids
#==============================================================================
#
#==============================================================================


def get_nearest_dwd_station(first_stn_id, coords_df_file,
                            x_col_name, y_col_name, index_col,
                            sep_type, neighbor_to_chose):
    ''' Find for one station, the closest neibhouring station'''
    # read df coordinates and get station ids, and x, y values
    in_coords_df, x_vals, y_vals, stn_ids = get_dwd_stns_coords(
        coords_df_file, x_col_name, y_col_name, index_col, sep_type)
    # make a tupples from the coordinates
    coords_tuples = np.array([(x, y) for x, y in zip(x_vals, y_vals)])
    # create a tree from coordinates
    points_tree = spatial.cKDTree(coords_tuples)

    xstn = in_coords_df.loc[int(first_stn_id), x_col_name]
    ystn = in_coords_df.loc[int(first_stn_id), y_col_name]

    distances, indices = points_tree.query([xstn, ystn], k=2)
    coords_nearest_nbr = coords_tuples[indices[neighbor_to_chose]]
    stn_near = str(stn_ids[indices[neighbor_to_chose]])
    distance_near = distances[neighbor_to_chose]

    return coords_nearest_nbr, stn_near, distance_near


#==============================================================================
#
#==============================================================================


def plt_bar_plot_2_stns(stn1_id, stn2_id, seperate_distance,
                        df1, df2, temp_freq, out_dir):
    ''' plot line plots between two stations'''
    print('plotting line plots')

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    time_vals = df1.index.to_pydatetime()
    time_arr = md.date2num(time_vals)

    ax.plot(time_arr, df1.values, c='darkblue',
            marker='o', markersize=2,
            alpha=0.25, label=stn1_id)
    ax.plot(time_arr, df2.values, c='darkred',
            marker='+', markersize=2,
            alpha=0.25, label=stn2_id)

    xfmt = md.DateFormatter('%Y-%m-%d')

    ax.xaxis.set_major_locator(MultipleLocator(35))

    ax.yaxis.set_major_locator(MultipleLocator(2))

    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(0.0, np.max([df1.values.max(), df2.values.max()]) + 1)
#     ax.set_xlabel('Time')
    ax.set_ylabel('Ppt in mm/%s' % temp_freq)
    ax.set_title("Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; " % (stn1_id, stn2_id,
                                      seperate_distance,
                                      temp_freq))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.xticks(rotation=45)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_lineplot_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Line plot')

    return
#==============================================================================
#
#==============================================================================


def plt_scatter_plot_2_stns(stn1_id, stn2_id, seperate_distance,
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
    _min = min(0, min(values_x.min(), values_y.min()))
    _max = max(values_x.max(), values_y.max())
    ax.plot([_min, _max], [_min, _max], c='k', linestyle='--', alpha=0.4)

    # set plot limit
    ax.set_xlim(-0.1, _max + 0.1)
    ax.set_ylim(-0.1, _max + 0.1)

    ax.xaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax.set_xlabel('Observed %s Rainfall in mm >= %0.1fmm Station Id: %s'
                  % (temp_freq, ppt_min_thr, stn1_id))
    ax.set_ylabel('Observed %s Rainfall in mm >= %0.1fmm Station Id: %s'
                  % (temp_freq, ppt_min_thr, stn2_id))
    ax.set_title("Stn: %s vs Stn: %s; \n Distance: %0.1f m; "
                 "Time Freq: %s; \n Pearson Cor=%0.3f; "
                 "Spearman Cor=%0.3f" % (stn1_id, stn2_id,
                                         seperate_distance,
                                         temp_freq, corr, rho))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_scatter_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Scatter')

    return

#==============================================================================
#
#==============================================================================


def plot_end_tail_cdf_2_stns(stn1_id, stn2_id, seperate_distance,
                             df1, df2, temp_freq, ppt_thr, out_dir):
    print('plotting CDF plots')
    values_x = df1.values.ravel()
    values_y = df2.values.ravel()

    xvals1, yvals1 = get_cdf_part_abv_thr(values_x, ppt_thr)
    xvals2, yvals2 = get_cdf_part_abv_thr(values_y, ppt_thr)

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)
    ax.plot(xvals1, yvals1, c='darkblue', marker='+', markersize=3,
            alpha=0.5, linestyle='--', label=stn1_id)
    ax.plot(xvals2, yvals2, c='darkred', marker='*', markersize=3,
            alpha=0.5, linestyle='dotted', label=stn2_id)
    ax.set_xlim(min(xvals1.min(), xvals2.min()) - 0.1,
                max(xvals1.max(), xvals2.max()) + 0.1)
    ax.set_ylim(min(yvals1.min(), yvals2.min()), 1)

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

    ax.set_xlabel('Observed %s Rainfall in mm stations Id: %s and %s'
                  % (temp_freq, stn1_id, stn2_id))
    ax.set_ylabel('Observed %s CDF stations Id: %s and %s'
                  % (temp_freq, stn1_id, stn2_id))
    ax.legend(loc='lower right')
    ax.set_title("CDF of extremes; Stn: %s vs Stn: %s;\n Distance: %0.1f m; \n"
                 "Time Freq: %s; Ppt thr: %.1f" % (stn1_id, stn2_id,
                                                   seperate_distance,
                                                   temp_freq, ppt_thr))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_cdf_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure CDF')

    return

#==============================================================================
#
#==============================================================================


def plot_normalized_ranked_stns(stn1_id, stn2_id, seperate_distance,
                                df1, df2, temp_freq, out_dir):
    print('plotting normalized Ranked plots')
    if isinstance(df1, pd.DataFrame):
        sorted_ranked_df1 = df1.rank(
            method='dense', na_option='keep')  # .sort_values(by=stn1_id)
    else:
        sorted_ranked_df1 = df1.rank(method='dense',
                                     na_option='keep')  # .sort_values()

    if isinstance(df2, pd.DataFrame):
        sorted_ranked_df2 = df2.rank(
            method='dense', na_option='keep')  # .sort_values(by=stn2_id)
    else:
        sorted_ranked_df2 = df2.rank(method='dense',
                                     na_option='keep')  # .sort_values()

    values_x = sorted_ranked_df1.values.ravel() / sorted_ranked_df1.values.max()
    values_y = sorted_ranked_df2.values.ravel() / sorted_ranked_df2.values.max()

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    ax.scatter(values_x, values_y, c='darkred', marker='.', alpha=0.35)

    # plot 45 deg line
    _min = min(values_x.min(), values_y.min())
    _max = max(values_x.max(), values_y.max())
    ax.plot([_min, _max], [_min, _max], c='k', linestyle='--', alpha=0.4)

    # set plot limit
    ax.set_xlim(-0.01, _max + 0.01)
    ax.set_ylim(-0.01, _max + 0.01)

    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_xlabel('Normalized Observed %s Ranks station Id: %s'
                  % (temp_freq, stn1_id))
    ax.set_ylabel('Normalized Observed %s Ranks station Id: %s'
                  % (temp_freq, stn2_id))
    ax.set_title("Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; "
                 % (stn1_id, stn2_id,
                     seperate_distance,
                     temp_freq))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_ranks_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure normalized Ranks')

    return
#==============================================================================
#
#==============================================================================


def plot_normalized_sorted_ranked_stns(stn1_id, stn2_id, seperate_distance,
                                       df1, df2, temp_freq, out_dir):
    print('plotting sorted normalized Ranked plots')
    if isinstance(df1, pd.DataFrame):
        sorted_ranked_df1 = df1.rank(
            method='dense', na_option='keep').sort_values(by=stn1_id)
    else:
        sorted_ranked_df1 = df1.rank(method='dense',
                                     na_option='keep').sort_values()

    if isinstance(df2, pd.DataFrame):
        sorted_ranked_df2 = df2.rank(
            method='dense', na_option='keep').sort_values(by=stn2_id)
    else:
        sorted_ranked_df2 = df2.rank(method='dense',
                                     na_option='keep').sort_values()

    values_x = sorted_ranked_df1.values.ravel() / sorted_ranked_df1.values.max()
    values_y = sorted_ranked_df2.values.ravel() / sorted_ranked_df2.values.max()

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    ax.scatter(values_x, values_y, c='darkred', marker='.', alpha=0.35)

    # plot 45 deg line
    _min = min(values_x.min(), values_y.min())
    _max = max(values_x.max(), values_y.max())
    ax.plot([_min, _max], [_min, _max], c='k', linestyle='--', alpha=0.4)

    # set plot limit
    ax.set_xlim(-0.01, _max + 0.01)
    ax.set_ylim(-0.01, _max + 0.01)

    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_xlabel('Normalized Observed %s Ranks station Id: %s'
                  % (temp_freq, stn1_id))
    ax.set_ylabel('Normalized Observed %s Ranks station Id: %s'
                  % (temp_freq, stn2_id))
    ax.set_title("Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; "
                 % (stn1_id, stn2_id,
                     seperate_distance,
                     temp_freq))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_sorted_ranks_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure sorted normalized Ranks')

    return

#==============================================================================
#
#==============================================================================


def plot_sorted_stns_vals(stn1_id, stn2_id, seperate_distance,
                          df1, df2, temp_freq, out_dir):
    print('plotting sorted plots')
    if isinstance(df1, pd.DataFrame):
        sorted_df1 = df1.sort_values(by=stn1_id)
    else:
        sorted_df1 = df1.sort_values()

    if isinstance(df2, pd.DataFrame):
        sorted_df2 = df2.sort_values(by=stn2_id)
    else:
        sorted_df2 = df2.sort_values()

    values_x = sorted_df1.values.ravel()
    values_y = sorted_df2.values.ravel()

    # calculate correlations (pearson and spearman)
    corr = pears(values_x, values_y)[0]
    rho = spr(values_x, values_y)[0]

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    ax.scatter(values_x, values_y, c='b', marker='+', alpha=0.5)

    # plot 45 deg line
    _min = min(values_x.min(), values_y.min())
    _max = max(values_x.max(), values_y.max())
    ax.plot([_min, _max], [_min, _max], c='k', linestyle='--', alpha=0.4)

    # set plot limit
    ax.set_xlim(-0.01, _max + 0.1)
    ax.set_ylim(-0.01, _max + 0.1)

    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_xlabel('Observed %s sorted Ppt station Id: %s'
                  % (temp_freq, stn1_id))
    ax.set_ylabel('Observed %s sorted Ppt Ranks station Id: %s'
                  % (temp_freq, stn2_id))
    ax.set_title("Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; \n Pearson Cor=%0.3f; "
                 "Spearman Cor=%0.3f" % (stn1_id, stn2_id,
                                         seperate_distance,
                                         temp_freq, corr, rho))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_sorted_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure Sorted dfs')

    return

#==============================================================================
#
#==============================================================================


def plot_p0_as_a_sequence_two_stns(stn_id,
                                   stn_2_id,
                                   min_dist,
                                   ppt_thrs_list,
                                   df_stn1,
                                   df_stn2,
                                   aggregation_frequencies_lst,
                                   out_dir):
    '''
    Plot P0 =probability that rainfall below a threshold as a function
    of time aggregations for every station and it's closest neighour
    '''
    print('Plotting P0 as a sequence')

    plt.ioff()

    colors = ['r', 'b', 'g', 'orange', 'k']
    df_p01_stn1 = pd.DataFrame(columns=aggregation_frequencies_lst,
                               index=ppt_thrs_list)
    df_p01_stn2 = pd.DataFrame(columns=aggregation_frequencies_lst,
                               index=ppt_thrs_list)
    fig = plt.figure(figsize=(16, 12), dpi=200)
    ax = fig.add_subplot(111)

    for i, ppt_thr in enumerate(ppt_thrs_list):
        print('Ppt Threshold is', ppt_thr)

        for temp_freq in aggregation_frequencies_lst:
            print('Time freq is', temp_freq)
            df_common1, df_common2 = resample_intersect_2_dfs(df_stn1,
                                                              df_stn2,
                                                              temp_freq)
            if (df_common1.values.shape[0] > 10 and
                    df_common2.values.shape[0] > 10):

                try:
                    p01 = calculate_probab_ppt_below_thr(
                        df_common1.values, ppt_thr)
                    p02 = calculate_probab_ppt_below_thr(
                        df_common2.values, ppt_thr)

                    df_p01_stn1.loc[ppt_thr, temp_freq] = p01
                    df_p01_stn2.loc[ppt_thr, temp_freq] = p02
                except Exception as msg:
                    print('error while calculating P0', msg, temp_freq)
                    continue

                print('plotting for Ppt thr', ppt_thr,
                      'temp freq ', temp_freq)

                ax.plot(df_p01_stn1.columns,
                        df_p01_stn1.loc[ppt_thr, :],
                        c=colors[i],
                        marker='o', linestyle='--',
                        linewidth=1,
                        alpha=0.5,
                        markersize=3,
                        label=str(ppt_thr) + ' mm')

                ax.plot(df_p01_stn2.columns,
                        df_p01_stn2.loc[ppt_thr, :], c=colors[i],
                        linewidth=1,
                        marker='+', linestyle='dotted', alpha=0.5,
                        markersize=3)

            else:
                print('empty df')
                break

    if df_p01_stn1.shape[0] > 0:

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

        ax.set_ylabel('P0')

        ax.set_title('P0 as a sequence %s vs %s \n distance: %0.1f m'
                     % (stn_id, stn_2_id, min_dist))

        plt.tight_layout()
        plt.grid(alpha=0.15)
        plt.savefig(os.path.join(out_dir,
                                 'P0_as_a_sequence_%s_%s.png'
                                 % (stn_id, stn_2_id)))
    else:
        print('empty df not plotting P0')

        print('Done plotting P0 as a sequence')
    return df_p01_stn1, df_p01_stn2

#==============================================================================
#
#==============================================================================


def plot_contingency_tables_as_a_sequence_two_stns(stn_id,
                                                   stn_2_id,
                                                   min_dist,
                                                   ppt_thrs_list,
                                                   df_stn1,
                                                   df_stn2,
                                                   aggregation_frequencies_lst,
                                                   out_dir):
    '''
    Do contingenccy tables a function of time and thresholds for 2 stations
    '''
    print('Plotting contingency_table a sequence')

    plt.ioff()

    for _, ppt_thr in enumerate(ppt_thrs_list):
        print('Ppt Threshold is', ppt_thr)

        for temp_freq in aggregation_frequencies_lst:

            print('Time freq is', temp_freq)

            fig = plt.figure(figsize=(16, 12), dpi=200)
            ax = fig.add_subplot(111)
            ax.set_aspect(1)
            df_common1, df_common2 = resample_intersect_2_dfs(df_stn1,
                                                              df_stn2,
                                                              temp_freq)
            if (df_common1.values.shape[0] > 10 and
                    df_common2.values.shape[0] > 10):
                print('getting percentages for contingency table')
                try:
                    (df_both_below_thr,
                     df_first_abv_second_below_thr,
                     df_first_below_second_abv_thr,
                     df_both_abv_thr) = constrcut_contingency_table(stn_id,
                                                                    stn_2_id,
                                                                    df_common1,
                                                                    df_common2,
                                                                    ppt_thr,
                                                                    ppt_thr)

                    conf_arr = np.array([[df_both_below_thr,
                                          df_first_abv_second_below_thr],
                                         [df_first_below_second_abv_thr,
                                          df_both_abv_thr]])

                except Exception as msg:
                    print('error while getting contingency table',
                          msg, temp_freq)
                    continue
            else:
                print('empty df values')
                break

            if conf_arr.shape[0] > 0:
                print('plotting for Ppt thr, temp frequency is', temp_freq)

                _ = sn.heatmap(conf_arr, annot=True, vmin=0.0,
                               vmax=100.0, fmt='.2f')

                ax.set_title('contingency table \n Rainfall Threshold %0.1fmm '
                             '%s vs %s \n distance: %0.1f m, time freq: %s'
                             % (ppt_thr, stn_id,
                                stn_2_id, min_dist, temp_freq))
                ax.set_xlabel('Station %s' % stn_2_id)
                ax.set_ylabel('Station %s' % stn_id)
                plt.savefig(
                    os.path.join(
                        out_dir,
                        'contingency_table_as_a_sequence_%s_%s_%s_%0.1fmm_.png'
                        % (stn_id, stn_2_id, temp_freq, ppt_thr)))

                plt.close(fig)
                print('Done plotting contingency_table as a sequence')

    return


#==============================================================================
#
#==============================================================================
def compare_two_dwd_stns(stns_ids):

    for iid in stns_ids:
        print('First Stn Id is', iid)
        try:
            idf1 = HDF52.get_pandas_dataframe(ids=[iid])
            idf1.dropna(axis=0, inplace=True)

            _, stn_near, distance_near = get_nearest_dwd_station(
                iid, coords_df_file, x_col_name, y_col_name, 3, ';',
                neighbor_to_chose)

            assert iid != stn_near, 'wrong neighbour selected'
            idf2 = HDF52.get_pandas_dataframe(ids=stn_near)
            print('Second Stn Id is', stn_near)

            out_save_dir = os.path.join(out_save_dir_orig,
                                        '%s_%s' % (iid, stn_near))

            if not os.path.exists(out_save_dir):
                os.mkdir(out_save_dir)

            for tem_freq in aggregation_frequencies:
                print('Aggregation is: ', tem_freq)
                df_common1, df_common2 = resample_intersect_2_dfs(idf1,
                                                                  idf2,
                                                                  tem_freq)

                if (df_common1.values.shape[0] > 1000 and
                        df_common2.values.shape[0] > 1000):
                    df_common1 = pd.DataFrame(
                        data=df_common1.values,
                        index=df_common1.index,
                        columns=[iid])

                    df_common2 = pd.DataFrame(
                        data=df_common2.values,
                        index=df_common2.index,
                        columns=[stn_near])
                    print('enough data are available for plotting')

                    try:
                        #====================================================

                        plt_bar_plot_2_stns(iid, stn_near, distance_near,
                                            df_common1, df_common2, tem_freq,
                                            out_save_dir)

                        plt_scatter_plot_2_stns(iid, stn_near, distance_near,
                                                df_common1, df_common2, ppt_thr_min,
                                                tem_freq,
                                                out_save_dir)
                        plot_end_tail_cdf_2_stns(iid, stn_near, distance_near,
                                                 df_common1, df_common2,
                                                 tem_freq, ppt_thr_min,
                                                 out_save_dir)
                        plot_normalized_ranked_stns(iid, stn_near,
                                                    distance_near,
                                                    df_common1, df_common2,
                                                    tem_freq,
                                                    out_save_dir)
                        plot_normalized_sorted_ranked_stns(iid, stn_near,
                                                           distance_near,
                                                           df_common1,
                                                           df_common2,
                                                           tem_freq,
                                                           out_save_dir)
                        plot_sorted_stns_vals(iid, stn_near, distance_near,
                                              df_common1, df_common2,
                                              tem_freq,
                                              out_save_dir)
                    except Exception as msg:
                        print('error while plotting', msg, tem_freq)
                        continue

                else:
                    print('Station is near but dont have enough data')
                    print('deleting out directory')
                    shutil.rmtree(out_save_dir, ignore_errors=True)
                    break

            if os.path.exists(out_save_dir):
                if (idf1.values.shape[0] > 1000 and
                        idf1.values.shape[0] > 1000):
                    # TODO: check is it works
                    ppt_thr = ppt_thrs_list[0]
                    print('Testing for Ppt Threshold of', ppt_thr)

                    for temp_freq in aggregation_frequencies:
                        print('Time freq is', temp_freq)
                        df_common1, df_common2 = resample_intersect_2_dfs(idf1,
                                                                          idf2,
                                                                          temp_freq)
                        if (df_common1.values.shape[0] > 10 and
                                df_common2.values.shape[0] > 10):
                            print(True, 'Plotting P0 and Contingency Tables')

                            # CHECK again with out dir
#                             plot_p0_as_a_sequence_two_stns(
#                                 iid,
#                                 stn_near,
#                                 distance_near,
#                                 ppt_thrs_list,
#                                 idf1,
#                                 idf2,
#                                 aggregation_frequencies,
#                                 out_save_dir)

                            plot_contingency_tables_as_a_sequence_two_stns(
                                iid,
                                stn_near,
                                distance_near,
                                ppt_thrs_list,
                                idf1,
                                idf2,
                                aggregation_frequencies,
                                out_save_dir)

                        else:
                            print('not enough data for P0 and Contingency Table')
                            break
                else:
                    print('not enough data for P0, Contnigency tables')
            if not os.listdir(out_save_dir):
                os.rmdir(out_save_dir)

        except Exception as msg:
            print(msg)

        break


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    HDF52 = HDF5(infile=path_to_ppt_hdf_data)
    ids = HDF52.get_all_ids()
    compare_two_dwd_stns(ids)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
