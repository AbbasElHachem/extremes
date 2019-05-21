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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as md

import scipy.spatial as spatial

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from pandas.plotting import register_matplotlib_converters

from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from _00_additional_functions import (resampleDf,
                                      get_cdf_part_abv_thr)

from b_get_data import HDF5

register_matplotlib_converters()

#==============================================================================
#
#==============================================================================
rc('font', size=13)
rc('font', family='serif')
rc('axes', labelsize=13)
rcParams['axes.labelpad'] = 13
majorLocator = MultipleLocator(1)

minorLocator = MultipleLocator(1)

path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')

coords_df_file = (r'X:\exchange\ElHachem\niederschlag_deutschland'
                  r'\1993_2016_5min_merge_nan.csv')
assert os.path.exists(coords_df_file), 'wrong DWD coords file'

out_save_dir = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_DWD_DWD')


if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)

x_col_name = 'Rechtswert'
y_col_name = 'Hochwert'

# threshold for CDF, consider only above thr, below is P0
ppt_thr = 1.

# till 1 day
aggregation_frequencies = ['5min', '10min', '15min', '30min', '60min', '90min',
                           '120min', '180min', '240min',  '360min',
                           '480min', '720min', '1440min']
# TODO: Fix for '300min', 520min

#==============================================================================
#
#==============================================================================


def get_dwd_stns_coords(coords_df_file, x_col_name, y_col_name):
    '''function used to return to coordinates dataframe of the DWD stations'''
    in_coords_df = pd.read_csv(coords_df_file, sep=';',
                               index_col=3, engine='c')
    stn_ids = in_coords_df.index
    x_vals = in_coords_df[x_col_name].values.ravel()
    y_vals = in_coords_df[y_col_name].values.ravel()
    return in_coords_df, x_vals, y_vals, stn_ids
#==============================================================================
#
#==============================================================================


def get_nearest_dwd_station(first_stn_id, coords_df_file,
                            x_col_name, y_col_name):
    ''' Find for one station, the closest neibhouring station'''
    # read df coordinates and get station ids, and x, y values
    in_coords_df, x_vals, y_vals, stn_ids = get_dwd_stns_coords(
        coords_df_file, x_col_name, y_col_name)
    # make a tupples from the coordinates
    coords_tuples = np.array([(x, y) for x, y in zip(x_vals, y_vals)])
    # create a tree from coordinates
    points_tree = spatial.cKDTree(coords_tuples)

    xstn = in_coords_df.loc[int(first_stn_id), x_col_name]
    ystn = in_coords_df.loc[int(first_stn_id), y_col_name]

    distances, indices = points_tree.query([xstn, ystn], k=2)
    coords_nearest_nbr = coords_tuples[indices[1]]
    stn_near = str(stn_ids[indices[1]])
    distance_near = distances[1]

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

    ax.plot(time_arr, df1.values, c='darkblue', marker='o', markersize=2,
            alpha=0.25, label=stn1_id)
    ax.plot(time_arr, df2.values, c='darkred', marker='+', markersize=2,
            alpha=0.25, label=stn2_id)

    xfmt = md.DateFormatter('%Y-%m-%d')

    ax.xaxis.set_major_locator(MultipleLocator(15))
    ax.yaxis.set_major_locator(MultipleLocator(2))

    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_ylim(0.0, np.max([df1.values.max(), df2.values.max()]) + 1)
    ax.set_xlabel('Observed %s Rainfall station Id: %s' % (temp_freq, stn1_id))
    ax.set_ylabel('Observed %s Rainfall station Id: %s' % (temp_freq, stn2_id))
    ax.set_title("Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; " % (stn1_id, stn2_id,
                                      seperate_distance,
                                      temp_freq))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.xticks(rotation=30)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_lineplot_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.close()
    print('Done saving figure Line plot')

    return
#==============================================================================
#
#==============================================================================


def plt_scatter_plot_2_stns(stn1_id, stn2_id, seperate_distance,
                            df1, df2, temp_freq, out_dir):
    ''' plot scatter plots between two stations'''
    print('plotting scatter plots')
    values_x = df1.values.ravel()
    values_y = df2.values.ravel()

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
    ax.set_xlim(-0.1, _max + 0.1)
    ax.set_ylim(-0.1, _max + 0.1)

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax.set_xlabel('Observed %s Rainfall in mm station Id: %s'
                  % (temp_freq, stn1_id))
    ax.set_ylabel('Observed %s Rainfall in mm station Id: %s'
                  % (temp_freq, stn2_id))
    ax.set_title("Stn: %s vs Stn: %s; \n Distance: %0.1f m; "
                 "Time Freq: %s; \n Pearson Cor=%0.3f; "
                 "Spearman Cor=%0.3f" % (stn1_id, stn2_id,
                                         seperate_distance,
                                         temp_freq, corr, rho))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_scatter_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.close()
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
    plt.close()
    print('Done saving figure CDF')

    return

#==============================================================================
#
#==============================================================================


def plot_ranked_stns(stn1_id, stn2_id, seperate_distance,
                     df1, df2, temp_freq, out_dir):
    print('plotting Ranked plots')
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

    values_x = sorted_ranked_df1.values.ravel()
    values_y = sorted_ranked_df2.values.ravel()

    # calculate correlations (pearson and spearman)
    corr = pears(values_x, values_y)[0]
    rho = spr(values_x, values_y)[0]

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)

    ax.scatter(values_x, values_y, c='darkred', marker='.', alpha=0.35)

    # plot 45 deg line
    _min = min(values_x.min(), values_y.min())
    _max = max(values_x.max(), values_y.max())
    ax.plot([_min, _max], [_min, _max], c='k', linestyle='--', alpha=0.4)

    # set plot limit
    ax.set_xlim(-0.1, _max + 0.1)
    ax.set_ylim(-0.1, _max + 0.1)

    ax.xaxis.set_major_locator(MultipleLocator(25))
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.set_xlabel('Observed %s Ranks station Id: %s' % (temp_freq, stn1_id))
    ax.set_ylabel('Observed %s Ranks station Id: %s' % (temp_freq, stn2_id))
    ax.set_title("Stn: %s vs Stn: %s;\n Distance: %0.1f m; "
                 "Time Freq: %s; \n Pearson Cor=%0.3f; "
                 "Spearman Cor=%0.3f" % (stn1_id, stn2_id,
                                         seperate_distance,
                                         temp_freq, corr, rho))

    ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.25)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '%s_ranks_stn_%s_vs_stn_%s_.png'
                             % (temp_freq, stn1_id, stn2_id)))
    plt.close()
    print('Done saving figure Ranks')

    return

#==============================================================================
#
#==============================================================================


def compare_cdf_two_dwd_stns(stns_ids):
    for iid in stns_ids[10:]:
        print('First Stn Id is', iid)
        try:
            idf1 = HDF52.get_pandas_dataframe(ids=[iid])
            idf1.dropna(axis=0, inplace=True)

            _, stn_near, distance_near = get_nearest_dwd_station(
                iid, coords_df_file, x_col_name, y_col_name)
            assert iid != stn_near, 'wrong neighbour selected'
            idf2 = HDF52.get_pandas_dataframe(ids=stn_near)
            print('Second Stn Id is', stn_near)
            for tem_freq in aggregation_frequencies:
                print('Aggregation is: ', tem_freq)
                df_resample1 = resampleDf(data_frame=idf1,
                                          temp_freq=tem_freq)
                df_resample2 = resampleDf(data_frame=idf2,
                                          temp_freq=tem_freq)
#                 df2 = select_df_within_period(df_resample2,
#                                               df_resample1.index[0],
#                                               df_resample1.index[-1])
#                 df1 = select_df_within_period(df_resample1,
#                                               df2.index[0],
#                                               df2.index[-1])
                idx_common = df_resample1.index.intersection(
                    df_resample2.index)
                df_common1 = df_resample1.loc[idx_common, :]
                df_common2 = df_resample2.loc[idx_common, :]
                if (df_common1.values.shape[0] > 0 and
                        df_common2.values.shape[0] > 0):
                    try:
                        plt_bar_plot_2_stns(iid, stn_near, distance_near,
                                            df_common1, df_common2, tem_freq,
                                            out_save_dir)

                        plt_scatter_plot_2_stns(iid, stn_near, distance_near,
                                                df_common1, df_common2,
                                                tem_freq,
                                                out_save_dir)
                        plot_end_tail_cdf_2_stns(iid, stn_near, distance_near,
                                                 df_common1, df_common2,
                                                 tem_freq, ppt_thr,
                                                 out_save_dir)
                        plot_ranked_stns(iid, stn_near, distance_near,
                                         df_common1, df_common2,
                                         tem_freq,
                                         out_save_dir)
                    except Exception as msg:
                        print('error while plotting', msg, tem_freq)
                        continue
                        # break
                else:
                    print('empty df')
                    continue
        except Exception as msg:
            print(msg)

        break


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    HDF52 = HDF5(infile=path_to_ppt_hdf_data)
    ids = HDF52.get_all_ids()
    compare_cdf_two_dwd_stns(ids)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
