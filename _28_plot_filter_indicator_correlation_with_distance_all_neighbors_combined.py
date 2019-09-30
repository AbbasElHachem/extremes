# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Plot Spearman Indicator Correlation with Distance
Purpose: Filter data based on spearman indicator correlation

Created on: 2019-08-10

After calculating for every Netatmo station the Spearman indicator rank 
correlation with either neigboring DWD Netatmo station (neighbors 1 to 5)
plot the result, x axis for distance and y axis for correlation
find the distribution of correaltion with increasing distance and
choice of neighbor.

This is used as a filter to remove the 'bad' Netamo staions,by fitting
for all nieghbor combination (distance, correlation) one function


Parameters
----------

Input Files
    DF file for each neighbor, index is station name, columns contain 
    distance to neighbor and rank spearman correlation
    
    Either Netatmo-DWD or Netatmo-Netatmo or DWD-DWD 
    for different percentage threshold and different time aggregations
Returns
-------

    Plot correlation with distance for different neighbors
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

from _00_additional_functions import (func, fit_curve_get_vals_below_curve,
                                      remove_all_low_corr_short_dist,
                                      read_filter_df_corr_return_stns_x_y_vals,
                                      gen_path_df_file)

plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 12})
#==============================================================================
#
#==============================================================================
main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_DWD_ppt_correlation_'
data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'

netatmo_path_acc = r'year_allyears_df_comparing_correlations_max_sep_dist_1000000_'
dwd_path_Acc = r'year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '90'
time_freq = '60min'

data_source0 = 'Netatmo'  # reference station 'Netatmo'
data_source = 'dwd'  # compare to station 'netatmo'

#==============================================================================
#
#==============================================================================
if data_source0 == 'Netatmo':
    df0 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 0)
    df1 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 1)
    df2 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 2)
    df3 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 3)
    df4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 4)
    df5 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 5)
    df6 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 6)
    df7 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 7)
    save_dir = data_dir_Netamto_dfs
if data_source0 == 'DWD':
    # for DWD stations neighbors start from 1 not 0 !
    df0 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 1)
    df1 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 2)
    df2 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 3)
    df3 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 4)
    df4 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 5)
    df5 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 6)
    df6 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 7)
    df7 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 8)
    save_dir = data_dir_DWD_dfs

#==============================================================================


#==============================================================================


s0, x0, y0, in_df0 = read_filter_df_corr_return_stns_x_y_vals(df0)
s1, x1, y1, in_df1 = read_filter_df_corr_return_stns_x_y_vals(df1)
s2, x2, y2, in_df2 = read_filter_df_corr_return_stns_x_y_vals(df2)
s3, x3, y3, in_df3 = read_filter_df_corr_return_stns_x_y_vals(df3)
s4, x4, y4, in_df4 = read_filter_df_corr_return_stns_x_y_vals(df4)
s5, x5, y5, in_df5 = read_filter_df_corr_return_stns_x_y_vals(df5)
s6, x6, y6, in_df6 = read_filter_df_corr_return_stns_x_y_vals(df6)
s7, x7, y7, in_df7 = read_filter_df_corr_return_stns_x_y_vals(df7)
#==============================================================================

# x0, y0, s0 = remove_all_low_corr_short_dist(x0, y0, 5e3, 0.6, s0)
x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
    x0, y0, 5e3, 0.5, s0, in_df0)

x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
    x1, y1, 5e3, 0.5, s1, in_df1)

x2_gd_corr, y2_gd_corr, s2_gd_corr, df2_gd_corr = remove_all_low_corr_short_dist(
    x2, y2, 5e3, 0.5, s2, in_df2)

x3_gd_corr, y3_gd_corr, s3_gd_corr, df3_gd_corr = remove_all_low_corr_short_dist(
    x3, y3, 5e3, 0.5, s3, in_df3)

x4_gd_corr, y4_gd_corr, s4_gd_corr, df4_gd_corr = remove_all_low_corr_short_dist(
    x4, y4, 5e3, 0.5, s4, in_df4)

x5_gd_corr, y5_gd_corr, s5_gd_corr, df5_gd_corr = remove_all_low_corr_short_dist(
    x5, y5, 5e3, 0.5, s5, in_df5)
x6_gd_corr, y6_gd_corr, s6_gd_corr, df6_gd_corr = remove_all_low_corr_short_dist(
    x6, y6, 5e3, 0.5, s6, in_df6)
x7_gd_corr, y7_gd_corr, s7_gd_corr, df7_gd_corr = remove_all_low_corr_short_dist(
    x7, y7, 5e3, 0.5, s7, in_df7)

# =============================================================================

# combined all stns, distances and correlation data into arrays
stns = np.concatenate((s0_gd_corr, s1_gd_corr, s2_gd_corr,
                       s3_gd_corr, s4_gd_corr, s5_gd_corr,
                       s6_gd_corr, s7_gd_corr))
xs = np.concatenate((x0_gd_corr, x1_gd_corr, x2_gd_corr,
                     x3_gd_corr, x4_gd_corr,
                     x5_gd_corr, x6_gd_corr, x7_gd_corr))
ys = np.concatenate((y0_gd_corr, y1_gd_corr, y2_gd_corr,
                     y3_gd_corr, y4_gd_corr,
                     y5_gd_corr, y6_gd_corr, y7_gd_corr))

dfs = pd.DataFrame(index=stns)
dfs['Distance to neighbor'] = xs
dfs['Bool_Spearman_Correlation'] = ys

# xs, ys, stns, dfs = remove_all_low_corr_short_dist(
#     xs, ys, 5e3, 0.5, stns, dfs)

# remove nans
stns = stns[np.where(ys >= 0)]
xs = xs[np.where(ys >= 0)]
ys = ys[ys >= 0]

# fit function to all data combined
(yfs, xvals_below_curves, yvals_below_curves,
 xvals_above_curves, yvals_above_curves,
 stns_below_curves, stns_above_curves) = fit_curve_get_vals_below_curve(
    xs, ys, func, stns, shift_per=0.1)

# =============================================================================

stns_to_keep = np.unique(stns_above_curves)
stns_to_remove = np.unique(stns_below_curves)

# TODO: NEED FIXING ASK PROF
# stns_keep_all = np.intersect1d(
#     np.intersect1d(np.intersect1d(np.intersect1d(np.intersect1d(
#         stns_to_keep, s0),
#         s1), s2), s3), s4)

# y0_gd_final = y0_gd_corr

stns_keep_all = np.intersect1d(stns_to_keep, s0_gd_corr)

# s0_gd_final = stns_keep_all[y0_gd_corr >= y0_gd_final]

stns_keep_al_sr = pd.Series(stns_keep_all)
stns_keep_al_sr.to_csv(
    (main_dir /
        r'plots_NetAtmo_ppt_DWD_ppt_correlation_' /
        (r'keep_stns_good_any_neighbor_%s_per_.csv' % percent)),
    sep=';')
# =============================================================================

x0_abv = df0_gd_corr.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y0_abv = df0_gd_corr.loc[stns_keep_all,
                         'Bool_Spearman_Correlation'].dropna().values

x1_abv = df1_gd_corr.loc[stns_keep_all,
                         'Distance to neighbor']  # .dropna().values
y1_abv = df1_gd_corr.loc[stns_keep_all,
                         'Bool_Spearman_Correlation']  # .dropna().values

x2_abv = df2_gd_corr.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y2_abv = df2_gd_corr.loc[stns_keep_all,
                         'Bool_Spearman_Correlation'].dropna().values

x3_abv = df3_gd_corr.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y3_abv = df3_gd_corr.loc[stns_keep_all,
                         'Bool_Spearman_Correlation'].dropna().values

x4_abv = df4_gd_corr.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y4_abv = df4_gd_corr.loc[stns_keep_all,
                         'Bool_Spearman_Correlation'].dropna().values

x5_abv = df5_gd_corr.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y5_abv = df5_gd_corr.loc[stns_keep_all,
                         'Bool_Spearman_Correlation'].dropna().values

x6_abv = df6_gd_corr.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y6_abv = df6_gd_corr.loc[stns_keep_all,
                         'Bool_Spearman_Correlation'].dropna().values

x7_abv = df7_gd_corr.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y7_abv = df7_gd_corr.loc[stns_keep_all,
                         'Bool_Spearman_Correlation'].dropna().values
# TODO: Add Additional filter on stations
# =============================================================================

plt.ioff()
plt.figure(figsize=(16, 12), dpi=300)

marker_size_abv_curve = 34
marker_size_below_curve = 28
marker_size_curve = 15
#
plt.scatter(x0_abv, y0_abv, c='r', alpha=0.5,
            marker='x', label=('First Neighbor Stn nbr %d / %d'
                               % (y0_abv.shape[0], y0_gd_corr.shape[0])),
            s=marker_size_abv_curve)
plt.scatter(x1_abv, y1_abv, c='b', alpha=0.5,
            marker='.', label=('Second Neighbor Stn nbr %d / %d'
                               % (y1_abv.shape[0], y1_gd_corr.shape[0])),
            s=marker_size_abv_curve)
plt.scatter(x2_abv, y2_abv, c='g', alpha=0.5,
            marker='d', label=('Third Neighbor Stn nbr %d / %d'
                               % (y2_abv.shape[0], y2_gd_corr.shape[0])),
            s=marker_size_abv_curve)
plt.scatter(x3_abv, y3_abv, c='darkorange', alpha=0.5,
            marker='1', label=('Fourth Neighbor Stn nbr %d / %d'
                               % (y3_abv.shape[0], y3_gd_corr.shape[0])),
            s=marker_size_abv_curve)
plt.scatter(x4_abv, y4_abv, c='m', alpha=0.5,
            marker='1', label=('Fifth Neighbor Stn nbr %d / %d'
                               % (y4_abv.shape[0], y4_gd_corr.shape[0])),
            s=marker_size_abv_curve)

plt.scatter(x5_abv, y5_abv, c='k', alpha=0.5,
            marker='1', label=('Sixfth Neighbor Stn nbr %d / %d'
                               % (y5_abv.shape[0], y5_gd_corr.shape[0])),
            s=marker_size_abv_curve)

plt.scatter(x6_abv, y6_abv, c='y', alpha=0.5,
            marker='1', label=('Seventh Neighbor Stn nbr %d / %d'
                               % (y6_abv.shape[0], y6_gd_corr.shape[0])),
            s=marker_size_abv_curve)

plt.scatter(x7_abv, y7_abv, c='c', alpha=0.5,
            marker='1', label=('Eigth Neighbor Stn nbr %d / %d'
                               % (y7_abv.shape[0], y7_gd_corr.shape[0])),
            s=marker_size_abv_curve)

# plt.scatter(xvals_above_curves, yvals_above_curves, s=marker_size_abv_curve,
#             label="Stns above Curve %d / %d" % (yvals_above_curves.shape[0],
#                                                 stns.shape[0]),
#             alpha=.25, c='b', marker='d')
#
# plt.scatter(xs, yfs, label="Polynomial function",
#             alpha=.25, c='m', marker='*', s=marker_size_curve)
#
#
# plt.scatter(xvals_below_curves, yvals_below_curves, s=marker_size_below_curve,
#             label="Below Curve %d / %d" % (yvals_below_curves.shape[0],
#                                            stns.shape[0]),
#             alpha=.25, c='k', marker='X')

#plt.scatter([],[], label='Number of stations %d' % dfs_all.shape[0], marker='.')
plt.xlim([0, x7.max() + 1000])
plt.xticks(np.arange(0, x7.max() + 1000, 5000))
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
plt.tight_layout()
# plt.axis('equal')
plt.title('Keeping %s %d / %d stations: Indicator correlation with distance'
          ' for upper %s percent of data values'
          % (data_source0, y0_abv.shape[0], s0.shape[0], percent))
plt.savefig(save_dir /
            (r'_%s_%s_%s_percent_indic_corr_freq_%s_all_neighbors_filtered.png'
             % (data_source0, data_source, percent, time_freq)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
# plt.show()
plt.close()
