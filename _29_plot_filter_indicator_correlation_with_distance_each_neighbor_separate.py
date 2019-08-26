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

from _00_additional_functions import (func, fit_curve_get_vals_below_curve,
                                      gen_path_df_file,
                                      read_filter_df_corr_return_stns_x_y_vals,
                                      remove_all_low_corr_short_dist)

plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})

#==============================================================================
#
#==============================================================================
main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_DWD_ppt_correlation_'
data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'

netatmo_path_acc = r'year_allyears_df_comparing_correlations_max_sep_dist_500000_'
dwd_path_Acc = r'year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '99'
time_freq = '60min'

data_source0 = 'Netatmo'  # reference station 'Netatmo'
data_source = 'dwd'  # compare to station 'netatmo'
#==============================================================================

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
    save_dir = data_dir_DWD_dfs

#==============================================================================


s0, x0, y0, in_df0 = read_filter_df_corr_return_stns_x_y_vals(df0)
s1, x1, y1, in_df1 = read_filter_df_corr_return_stns_x_y_vals(df1)
s2, x2, y2, in_df2 = read_filter_df_corr_return_stns_x_y_vals(df2)
s3, x3, y3, in_df3 = read_filter_df_corr_return_stns_x_y_vals(df3)
s4, x4, y4, in_df4 = read_filter_df_corr_return_stns_x_y_vals(df4)

#==============================================================================

# x0, y0, s0 = remove_all_low_corr_short_dist(x0, y0, 5e3, 0.6, s0)
x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
    x1, y1, 5e3, 0.6, s1, in_df1)

x2_gd_corr, y2_gd_corr, s2_gd_corr, df2_gd_corr = remove_all_low_corr_short_dist(
    x2, y2, 5e3, 0.6, s2, in_df2)

x3_gd_corr, y3_gd_corr, s3_gd_corr, df3_gd_corr = remove_all_low_corr_short_dist(
    x3, y3, 5e3, 0.6, s3, in_df3)

x4_gd_corr, y4_gd_corr, s4_gd_corr, df4_gd_corr = remove_all_low_corr_short_dist(
    x4, y4, 5e3, 0.6, s4, in_df4)


# =============================================================================

# %% apply filter for every neighbor seperatly

#x, y, func, stns, df_data

(y_fitted_shifted0, xvals_below_curve0,
 yvals_below_curve0, xvals_above_curve0,
 yvals_above_curve0, stnbelow0, stnabove0) = fit_curve_get_vals_below_curve(
    x=x0, y=y0, func=func, stns=s0)

(y_fitted_shifted1, xvals_below_curve1,
 yvals_below_curve1, xvals_above_curve1,
 yvals_above_curve1, stnbelow1, stnabove1) = fit_curve_get_vals_below_curve(
    x=x1_gd_corr, y=y1_gd_corr, func=func, stns=s1_gd_corr)

(y_fitted_shifted2, xvals_below_curve2,
 yvals_below_curve2, xvals_above_curve2,
 yvals_above_curve2, stnbelow2, stnabove2) = fit_curve_get_vals_below_curve(
    x=x2_gd_corr, y=y2_gd_corr, func=func, stns=s2_gd_corr)

(y_fitted_shifted3, xvals_below_curve3,
 yvals_below_curve3, xvals_above_curve3,
 yvals_above_curve3, stnbelow3, stnabove3) = fit_curve_get_vals_below_curve(
    x=x3_gd_corr, y=y3_gd_corr, func=func, stns=s3_gd_corr)

(y_fitted_shifted4, xvals_below_curve4,
 yvals_below_curve4, xvals_above_curve4,
 yvals_above_curve4, stnbelow4, stnabove4) = fit_curve_get_vals_below_curve(
    x=x4_gd_corr, y=y4_gd_corr, func=func, stns=s4_gd_corr)

# =============================================================================


# TODO: CHECK AGAIN intersect all stations
# https://www.python-course.eu/polynomial_class_in_python.php
stns_keep_all = np.intersect1d(np.intersect1d(np.intersect1d(np.intersect1d(
    stnabove0, stnabove1),
    stnabove2), stnabove3), stnabove4)

stns_keep_al_sr = pd.DataFrame(data=stnabove0,
                               columns=['Stations'])
stns_keep_al_sr.to_csv(
    (save_dir /
        (r'keep_stns_1st_neighbor_%s_per_%s_.csv'
         % (percent, time_freq))),
    sep=';')
# # %%
# x0_abv = in_df0.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
# y0_abv = in_df0.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values
#
# x1_abv = in_df1.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
# y1_abv = in_df1.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values
#
# x2_abv = in_df2.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
# y2_abv = in_df2.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values
#
# x3_abv = in_df3.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
# y3_abv = in_df3.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values
#
# x4_abv = in_df4.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
# y4_abv = in_df4.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

x0_abv = in_df0.loc[stnabove0, 'Distance to neighbor'].dropna().values
y0_abv = in_df0.loc[stnabove0, 'Bool_Spearman_Correlation'].dropna().values

x1_abv = in_df1.loc[stnabove1, 'Distance to neighbor'].dropna().values
y1_abv = in_df1.loc[stnabove1, 'Bool_Spearman_Correlation'].dropna().values

x2_abv = in_df2.loc[stnabove2, 'Distance to neighbor'].dropna().values
y2_abv = in_df2.loc[stnabove2, 'Bool_Spearman_Correlation'].dropna().values

x3_abv = in_df3.loc[stnabove3, 'Distance to neighbor'].dropna().values
y3_abv = in_df3.loc[stnabove3, 'Bool_Spearman_Correlation'].dropna().values

x4_abv = in_df4.loc[stnabove4, 'Distance to neighbor'].dropna().values
y4_abv = in_df4.loc[stnabove4, 'Bool_Spearman_Correlation'].dropna().values

#==============================================================================
# PLOT
#==============================================================================
plt.ioff()
# plt.figure(figsize=(36, 18), dpi=300)
#
# marker_size_abv_curve = 100
# marker_size_below_curve = 95
# marker_size_curve = 55

plt.figure(figsize=(16, 12), dpi=300)

marker_size_abv_curve = 34
marker_size_below_curve = 28
marker_size_curve = 15

plt.scatter(x0_abv, y0_abv, c='r', alpha=0.5,
            marker='x', label='First Neighbor Stn nbr %d' % y0_abv.shape[0],
            s=marker_size_abv_curve)
plt.scatter(x1_abv, y1_abv, c='b', alpha=0.5,
            marker='.', label='Second Neighbor Stn nbr %d' % y1_abv.shape[0],
            s=marker_size_abv_curve)
plt.scatter(x2_abv, y2_abv, c='g', alpha=0.5,
            marker='d', label='Third Neighbor Stn nbr %d' % y2_abv.shape[0],
            s=marker_size_abv_curve)
plt.scatter(x3_abv, y3_abv, c='darkorange', alpha=0.5,
            marker='*', label='Fourth NeighborStn nbr %d' % y3_abv.shape[0],
            s=marker_size_abv_curve)
plt.scatter(x4_abv, y4_abv, c='m', alpha=0.5,
            marker='+', label='Fifth Neighbor Stn nbr %d' % y4_abv.shape[0],
            s=marker_size_abv_curve)

plt.scatter(xvals_below_curve0, yvals_below_curve0, c='grey', alpha=0.65,
            marker='x', s=marker_size_below_curve)
plt.scatter(xvals_below_curve1, yvals_below_curve1, c='grey', alpha=0.65,
            marker='.',  s=marker_size_below_curve)
plt.scatter(xvals_below_curve2, yvals_below_curve2, c='grey', alpha=0.65,
            marker='d',  s=marker_size_below_curve)
plt.scatter(xvals_below_curve3, yvals_below_curve3, c='grey', alpha=0.65,
            marker='*',  s=marker_size_below_curve)
plt.scatter(xvals_below_curve4, yvals_below_curve4, c='grey', alpha=0.65,
            marker='+', s=marker_size_below_curve)

plt.scatter(x0, y_fitted_shifted0, c='darkred', alpha=0.25,
            marker='.',  s=marker_size_curve)  # label='Fitted curve 1',

plt.scatter(x1_gd_corr, y_fitted_shifted1, c='darkblue', alpha=0.25,
            marker='.', s=marker_size_curve)  # label='Fitted curve 2',

plt.scatter(x2_gd_corr, y_fitted_shifted2, c='darkgreen', alpha=0.25,
            marker='.', s=marker_size_curve)  # label='Fitted curve 3',

plt.scatter(x3_gd_corr, y_fitted_shifted3, c='darkorange', alpha=0.25,
            marker='.',  s=marker_size_curve)  # label='Fitted curve 4',

plt.scatter(x4_gd_corr, y_fitted_shifted4, c='purple', alpha=0.25,
            marker='.', s=marker_size_curve)  # label='Fitted curve 5',
# plt.show()

plt.xlim([0, max([x3.max(), x4.max()]) + 1000])
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)

plt.title('Keeping %s %d stations: Indicator correlation with distance'
          ' for upper %s percent for %s data values (each neighbor)'
          % (data_source0, stns_keep_al_sr.shape[0], percent, time_freq))
plt.savefig(save_dir /
            (r'_%s_%s_%s_percent_indic_corr_freq_%s_filtered_each_neighbors_data.png'
             % (data_source0, data_source, percent, time_freq)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
