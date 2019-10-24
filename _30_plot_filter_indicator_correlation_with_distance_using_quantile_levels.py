# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Plot Spearman Indicator Correlation with Distance
Purpose: Filter data based on spearman indicator correlation using each 
    neighbor seperatly

Created on: 2019-08-10

After calculating for every Netatmo station the Spearman indicator rank 
correlation with either neigboring DWD Netatmo station (neighbors 1 to 5)
plot the result, x axis for distance and y axis for correlation
find the distribution of correaltion with increasing distance and
choice of neighbor.

This is used as a filter to remove the 'bad' Netamo staions, by fitting
for each nieghbor combination (distance, correlation) a function


Parameters
----------

Input Files
    DF file for each neighbor, index is station name, columns contain 
    distance to neighbor and rank spearman correlation
    
    Either Netatmo-DWD or Netatmo-Netatmo or DWD-DWD 
    for different percentage threshold and different time aggregations
Returns
-------

    Plot correlation with distance for different neighbors and fitted curves
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path
from random import shuffle

from _00_additional_functions import (func, find_nearest,
                                      gen_path_df_file,
                                      read_filter_df_corr_return_stns_x_y_vals,
                                      remove_all_low_corr_short_dist,
                                      chunks)

plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})

#==============================================================================
#
#==============================================================================
main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_DWD_ppt_correlation_'
data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'

netatmo_path_acc = r'pearson_year_allyears_df_comparing_correlations_max_sep_dist_30000_'
dwd_path_Acc = r'year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '97_0'
time_freq = '60min'

data_source0 = 'Netatmo'  # reference station 'Netatmo'

data_source = 'dwd'  # compare to station 'netatmo'

remove_upper_limit = False

plot_dwd_on_top = True

method1 = False
method2 = True
#==============================================================================

#==============================================================================
if data_source0 == 'Netatmo':
    df0 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 0)
    df1 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 1)
    df2 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 2)
#     df3 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 3)
#     df4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 4)
#     df5 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 5)
#     df6 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 6)
#     df7 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 7)
    save_dir = data_dir_Netamto_dfs
if data_source0 == 'DWD':
    percent_dwd = percent.replace('_0', '')
    # for DWD stations neighbors start from 1 not 0 !
    df0 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent_dwd, 1)
    df1 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent_dwd, 2)
    df2 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent_dwd, 3)
#     df3 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 4)
#     df4 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 5)
#     df5 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 6)
#     df6 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 7)
#     df7 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 8)
    save_dir = data_dir_DWD_dfs

if plot_dwd_on_top:
    percent_dwd = percent.replace('_0', '')
    # for DWD stations neighbors start from 1 not 0 !
    df0_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent_dwd, 1)
    df1_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent_dwd, 2)
    df2_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent_dwd, 3)
#     df3 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 4)
#     df4 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 5)
#     df5 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 6)
#     df6 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 7)
#     df7 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 8)

#==============================================================================


s0, x0, y0, in_df0 = read_filter_df_corr_return_stns_x_y_vals(df0)
s1, x1, y1, in_df1 = read_filter_df_corr_return_stns_x_y_vals(df1)

# remove above 20km
x1 = x1[np.where(0 < x1) and np.where(x1 <= 20000)]
x1 = x1[x1 > 0]

s1 = s1[np.where(0 < x1) and np.where(x1 <= 20000)]
y1 = y1[np.where(0 < x1) and np.where(x1 <= 20000)]
in_df1 = in_df1.loc[in_df1.index.intersection(s1), :]

# s2, x2, y2, in_df2 = read_filter_df_corr_return_stns_x_y_vals(df2)
# s3, x3, y3, in_df3 = read_filter_df_corr_return_stns_x_y_vals(df3)
# s4, x4, y4, in_df4 = read_filter_df_corr_return_stns_x_y_vals(df4)
# s5, x5, y5, in_df5 = read_filter_df_corr_return_stns_x_y_vals(df5)
# s6, x6, y6, in_df6 = read_filter_df_corr_return_stns_x_y_vals(df6)
# s7, x7, y7, in_df7 = read_filter_df_corr_return_stns_x_y_vals(df7)

if plot_dwd_on_top:
    s0_dwd, x0_dwd, y0_dwd, in_df0_dwd = read_filter_df_corr_return_stns_x_y_vals(
        df0_dwd)
    s1_dwd, x1_dwd, y1_dwd, in_df1_dwd = read_filter_df_corr_return_stns_x_y_vals(
        df1_dwd)
#==============================================================================

# x0, y0, s0 = remove_all_low_corr_short_dist(x0, y0, 5e3, 0.6, s0)

if percent == '95_0':
    x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
        x0, y0, 5e3, 0.5, s0, in_df0)

    x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
        x1, y1, 5e3, 0.5, s1, in_df1)

if percent == '97_0':
    x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
        x0, y0, 5e3, 0.5, s0, in_df0)

    x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
        x1, y1, 5e3, 0.5, s1, in_df1)

if percent == '98_0':
    x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
        x0, y0, 5e3, 0.4, s0, in_df0)

    x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
        x1, y1, 5e3, 0.4, s1, in_df1)

if percent == '99_0':
    x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
        x0, y0, 5e3, 0.1, s0, in_df0)

    x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
        x1, y1, 5e3, 0.025, s1, in_df1)
# x2_gd_corr, y2_gd_corr, s2_gd_corr, df2_gd_corr = remove_all_low_corr_short_dist(
#     x2, y2, 5e3, 0.5, s2, in_df2)
#
# x3_gd_corr, y3_gd_corr, s3_gd_corr, df3_gd_corr = remove_all_low_corr_short_dist(
#     x3, y3, 5e3, 0.5, s3, in_df3)
#
# x4_gd_corr, y4_gd_corr, s4_gd_corr, df4_gd_corr = remove_all_low_corr_short_dist(
#     x4, y4, 5e3, 0.5, s4, in_df4)
#
# x5_gd_corr, y5_gd_corr, s5_gd_corr, df5_gd_corr = remove_all_low_corr_short_dist(
#     x5, y5, 5e3, 0.5, s5, in_df5)
# x6_gd_corr, y6_gd_corr, s6_gd_corr, df6_gd_corr = remove_all_low_corr_short_dist(
#     x6, y6, 5e3, 0.5, s6, in_df6)
# x7_gd_corr, y7_gd_corr, s7_gd_corr, df7_gd_corr = remove_all_low_corr_short_dist(
#     x7, y7, 5e3, 0.5, s7, in_df7)
# =============================================================================


#x, y, func, stns, df_data

# all_netatmo_dwd_stns_pairs = list(np.unique(
#     np.intersect1d(s1_gd_corr.values, s0_gd_corr.values)))
# print('Shuffling unique stns ids', len(all_netatmo_dwd_stns_pairs))
#
#
# shuffle(all_netatmo_dwd_stns_pairs)
# shuffled_netatmo_dwd_stns_grps = np.array(
#     list(chunks(all_netatmo_dwd_stns_pairs,
#                 int(np.round(len(all_netatmo_dwd_stns_pairs) / y0_dwd.size,
#                              0)))))

all_stns = list(np.hstack((s0_gd_corr.values, s1_gd_corr.values)))
# print(shuffled_netatmo_dwd_stns_grps.size, ' groups of unique stns ')

if method1:
    stns_keep = []

    corr_keep0, corr_keep1 = [], []
    dist_keep0, dist_keep1 = [], []

    for netatmo_stn_one in all_stns:
        # print(netatmo_stn_one)
        #     for netatmo_stn_one in netatmo_stns:

        try:

            corr_stn0 = df0_gd_corr.loc[netatmo_stn_one,
                                        'Bool_Spearman_Correlation']
            distance_stn0 = df0_gd_corr.loc[netatmo_stn_one,
                                            'Distance to neighbor']

            corr_stn1 = df1_gd_corr.loc[netatmo_stn_one,
                                        'Bool_Spearman_Correlation']
            distance_stn1 = df1_gd_corr.loc[netatmo_stn_one,
                                            'Distance to neighbor']

            x0_dwd_nearest = find_nearest(x0_dwd, distance_stn0)
            y0_dwd_nearest = y0_dwd[np.where(x0_dwd == x0_dwd_nearest)][0]

            x1_dwd_nearest = find_nearest(x1_dwd, distance_stn1)
            y1_dwd_nearest = y1_dwd[np.where(x1_dwd == x1_dwd_nearest)][0]

            if corr_stn0 < y0_dwd_nearest.min() or corr_stn1 < y1_dwd_nearest.min():
                print('Station is removed')

            if corr_stn0 >= y0_dwd_nearest.min() and corr_stn1 >= y1_dwd_nearest.min():

                if (netatmo_stn_one in s0_gd_corr) and (
                        netatmo_stn_one in s1_gd_corr):
                    dist_keep0.append(distance_stn0)
                    dist_keep1.append(distance_stn1)

                    corr_keep0.append(corr_stn0)
                    corr_keep1.append(corr_stn1)

                    if (netatmo_stn_one not in stns_keep):
                        print('adding stn', netatmo_stn_one)
                        stns_keep.append(netatmo_stn_one)

    #             corrs_stn = [all_corrs[i] for i, stn in enumerate(all_stns) if
    #                          stn == netatmo_stn_one]
    #
    #             dist_stn = [all_dists[i] for i, stn in enumerate(all_stns) if
    #                         stn == netatmo_stn_one]
    #
    #             idx_min_corr = np.argmin(corrs_stn)
    #             dist_min = dist_stn[idx_min_corr]
    #             if min(corrs_stn) >= dwd_min_corr:
    #
    #                 corr_keep.append(min(corrs_stn))
    #
    #                 dist_keep.append(dist_min)
        except Exception as msg:
            print(msg)
            continue
    #     pass
    print('done getting stations')

    print('keeping', len(stns_keep))

#==============================================================================
#
#==============================================================================

if method2:

    dwd_min_corr = min(y0_dwd.min(), y1_dwd.min())
    # %% apply filter for every neighbor seperatly

    all_dists = np.hstack((x0_gd_corr, x1_gd_corr))
    all_corrs = np.hstack((y0_gd_corr, y1_gd_corr))

    stns_keep = []

    corr_keep0, corr_keep1 = [], []
    dist_keep0, dist_keep1 = [], []

    for netatmo_stn_one in all_stns:
        # print(netatmo_stn_one)
        #     for netatmo_stn_one in netatmo_stns:

        try:

            corr_stn0 = df0_gd_corr.loc[netatmo_stn_one,
                                        'Bool_Spearman_Correlation']
            distance_stn0 = df0_gd_corr.loc[netatmo_stn_one,
                                            'Distance to neighbor']

            corr_stn1 = df1_gd_corr.loc[netatmo_stn_one,
                                        'Bool_Spearman_Correlation']
            distance_stn1 = df1_gd_corr.loc[netatmo_stn_one,
                                            'Distance to neighbor']

            x0_dwd_nearest = find_nearest(x0_dwd, distance_stn0)
            y0_dwd_nearest = y0_dwd[np.where(x0_dwd == x0_dwd_nearest)][0]

            x1_dwd_nearest = find_nearest(x1_dwd, distance_stn1)
            y1_dwd_nearest = y1_dwd[np.where(x1_dwd == x1_dwd_nearest)][0]

            if corr_stn0 < y0_dwd.min() or corr_stn1 < y1_dwd.min():
                print('Station is removed')

            if corr_stn0 >= y0_dwd.min() and corr_stn1 >= y1_dwd.min():

                if (netatmo_stn_one in s0_gd_corr) and (
                        netatmo_stn_one in s1_gd_corr):

                    if (netatmo_stn_one not in stns_keep):
                        print('adding stn', netatmo_stn_one)
                        stns_keep.append(netatmo_stn_one)
                        dist_keep0.append(distance_stn0)
                        dist_keep1.append(distance_stn1)

                        corr_keep0.append(corr_stn0)
                        corr_keep1.append(corr_stn1)

    #             corrs_stn = [all_corrs[i] for i, stn in enumerate(all_stns) if
    #                          stn == netatmo_stn_one]
    #
    #             dist_stn = [all_dists[i] for i, stn in enumerate(all_stns) if
    #                         stn == netatmo_stn_one]
    #
    #             idx_min_corr = np.argmin(corrs_stn)
    #             dist_min = dist_stn[idx_min_corr]
    #             if min(corrs_stn) >= dwd_min_corr:
    #
    #                 corr_keep.append(min(corrs_stn))
    #
    #                 dist_keep.append(dist_min)
        except Exception as msg:
            print(msg)
            continue

    print('keeping', len(stns_keep))
# x_stns_keep0 = in_df0.loc[stns_keep, 'Distance to neighbor'].values
# y_stns_keep0 = in_df0.loc[stns_keep, 'Bool_Spearman_Correlation'].values
#
# x_stns_keep1 = in_df1.loc[stns_keep, 'Distance to neighbor'].values
# y_stns_keep1 = in_df1.loc[stns_keep, 'Bool_Spearman_Correlation'].values

#==============================================================================
#
#==============================================================================

# # stns_keep_all_final = stn0_below_22
# stns_keep_al_sr = pd.DataFrame(data=stns_keep_all_final_new,
#                                columns=['Stations'])
#
# stns_keep_al_sr.to_csv(
#     (save_dir /
#         (r'keep_stns_all_neighbor_%s_per_%s_s0_1.csv'
#          % (percent, time_freq))),
#     sep=';')
#

#==============================================================================
#
#==============================================================================
plt.ioff()
# plt.figure(figsize=(36, 18), dpi=300)
#
# marker_size_abv_curve = 100
# marker_size_below_curve = 95
# marker_size_curve = 55

plt.figure(figsize=(16, 12), dpi=300)

marker_size_abv_curve = 32
marker_size_below_curve = 24
marker_size_curve = 12

# plt.scatter(dist_keep, corr_keep, c='r', alpha=0.5,
#             marker='x', label=('Remaining Stn nbr %d / %d'
#                                % (len(corr_keep), len(all_stns))),  # .shape[0])),
#             s=marker_size_abv_curve)
plt.scatter(dist_keep0, corr_keep0, c='r', alpha=0.5,
            marker='x', label=('Remaining Stn nbr %d / %d'
                               % (len(corr_keep0), len(s0_gd_corr))),
            s=marker_size_abv_curve)
plt.scatter(dist_keep1, corr_keep1, c='b', alpha=0.5,
            marker='.', label=('Remaining Stn nbr %d / %d'
                               % (len(corr_keep1), len(s1_gd_corr))),
            s=marker_size_abv_curve)
if plot_dwd_on_top:
    plt.scatter(x0_dwd, y0_dwd, c='k', alpha=0.5,
                marker='D', label=('DWD First Neighbor Stn nbr %d'
                                   % (y0_dwd.shape[0])),
                s=marker_size_abv_curve)
    plt.scatter(x1_dwd, y1_dwd, c='k', alpha=0.5,
                marker='X', label=('DWD Second Neighbor Stn nbr %d '
                                   % (y1_dwd.shape[0])),
                s=marker_size_abv_curve)
#
#
plt.xlim([0, max([x0_gd_corr.max(), x1_gd_corr.max()]) + 1000])
plt.xticks(np.arange(0, x1_gd_corr.max() + 1000, 5000))
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)

plt.title('Keeping %s %d stations: Indicator correlation with distance'
          ' for upper %s percent for %s data values'
          % (data_source0,  len(corr_keep0),
              percent, time_freq))
plt.savefig(save_dir /  # filtered_2_neighbors_data
            (r'_%s_%s_%s_percent_indic_corr_freq_%s_filtered_pearson2.png'
             % (data_source0, data_source, percent, time_freq)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
