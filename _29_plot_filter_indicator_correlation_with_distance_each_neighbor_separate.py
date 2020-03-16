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

plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'axes.labelsize': 26})

#==============================================================================
#
#==============================================================================
main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

# for BW
# data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_DWD_ppt_correlation_'
# data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'

# for RH
data_dir_Netamto_dfs = Path(
    r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\indicator_correlation')

data_dir_DWD_dfs = data_dir_Netamto_dfs
netatmo_path_acc_b4 = r'pearson_year_allyears_df_comparing_correlations_max_sep_dist_1000000_'

netatmo_path_acc = r'pearson_year_allyears_df_comparing_correlations_max_sep_dist_1000000_'

dwd_path_Acc = r'pearson_year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '99_5'
time_freq = '60min'

data_source0 = 'Netatmo'  # reference station 'Netatmo'

data_source = 'dwd'  # compare to station 'netatmo'

remove_upper_limit = False

plot_dwd_on_top = True

shift_by_percent = 0
#==============================================================================

#==============================================================================
# if data_source0 == 'Netatmo':
#     df0 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 0)
#     df1 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 1)
#     df2 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 2)
#
#     save_dir = data_dir_Netamto_dfs
#

# if data_source0 == 'DWD':
#     #     percent_dwd = percent.replace('_0', '')
#     # for DWD stations neighbors start from 1 not 0 !
#     df0 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 1)
#     df1 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 2)
#     df2 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 3)
#
#     save_dir = data_dir_DWD_dfs

if plot_dwd_on_top:
    # percent_dwd = percent.replace('_0', '')
    # for DWD stations neighbors start from 1 not 0 !
    df0_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent, 1)
    df1_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent, 2)
    df2_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent, 3)

df0_b4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_b4, time_freq,
                          data_source, percent, 0)
df1_b4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_b4, time_freq,
                          data_source, percent, 1)
df2_b4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_b4, time_freq,
                          data_source, percent, 2)

df0_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 1)
df1_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 2)
df2_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 3)

save_dir = data_dir_Netamto_dfs
#==============================================================================

in_df0_dwd = pd.read_csv(df0_dwd, index_col=0, sep=';').dropna(how='all')
in_df1_dwd = pd.read_csv(df1_dwd, index_col=0, sep=';').dropna(how='all')
in_df2_dwd = pd.read_csv(df2_dwd, index_col=0, sep=';').dropna(how='all')

s0, x0, y0, in_df0 = read_filter_df_corr_return_stns_x_y_vals(df0_b4)
# s1, x1, y1, in_df1 = read_filter_df_corr_return_stns_x_y_vals(df1)

# remove above 20km
# x1 = x1[np.where(0 < x1) and np.where(x1 <= 20000)]
# x1 = x1[x1 > 0]
#
# s1 = s1[np.where(0 < x1) and np.where(x1 <= 20000)]
# y1 = y1[np.where(0 < x1) and np.where(x1 <= 20000)]
# in_df1 = in_df1.loc[in_df1.index.intersection(s1), :]

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

# if percent == '95_0':
#     x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
#         x0, y0, 5e3, 0.5, s0, in_df0)

#     x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
#         x1, y1, 5e3, 0.5, s1, in_df1)

# if percent == '97_0':
#     x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
#         x0, y0, 5e3, 0.5, s0, in_df0)
#
#     x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
#         x1, y1, 5e3, 0.5, s1, in_df1)
#
# if percent == '98_0':
#     x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
#         x0, y0, 5e3, 0.4, s0, in_df0)
#
#     x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
#         x1, y1, 5e3, 0.4, s1, in_df1)
#
# if percent == '99_0':
#     x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
#         x0, y0, 5e3, 0.1, s0, in_df0)
#
#     x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
#         x1, y1, 5e3, 0.025, s1, in_df1)
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

# %% apply filter for every neighbor seperatly

#x, y, func, stns, df_data

(y_fitted_shifted0, xvals_below_curve0,
 yvals_below_curve0, xvals_above_curve0,
 yvals_above_curve0, stnbelow0, stnabove0) = fit_curve_get_vals_below_curve(
    x=x0, y=y0, func=func, stns=s0, shift_per=shift_by_percent / 100)

# (y_fitted_shifted1, xvals_below_curve1,
#  yvals_below_curve1, xvals_above_curve1,
#  yvals_above_curve1, stnbelow1, stnabove1) = fit_curve_get_vals_below_curve(
#     x=x1_gd_corr, y=y1_gd_corr, func=func, stns=s1_gd_corr, shift_per=0.001)

# if percent == '98_0':
#     (y_fitted_shifted0, xvals_below_curve0,
#      yvals_below_curve0, xvals_above_curve0,
#      yvals_above_curve0, stnbelow0, stnabove0) = fit_curve_get_vals_below_curve(
#         x=x0_gd_corr, y=y0_gd_corr, func=func, stns=s0_gd_corr, shift_per=0.15)
#
#     (y_fitted_shifted1, xvals_below_curve1,
#      yvals_below_curve1, xvals_above_curve1,
#      yvals_above_curve1, stnbelow1, stnabove1) = fit_curve_get_vals_below_curve(
# x=x1_gd_corr, y=y1_gd_corr, func=func, stns=s1_gd_corr, shift_per=0.05)

# if percent == '99_0':
#     (y_fitted_shifted0, xvals_below_curve0,
#      yvals_below_curve0, xvals_above_curve0,
#      yvals_above_curve0, stnbelow0, stnabove0) = fit_curve_get_vals_below_curve(
#         x=x0_gd_corr, y=y0_gd_corr, func=func, stns=s0_gd_corr, shift_per=0.05)
#
#     (y_fitted_shifted1, xvals_below_curve1,
#      yvals_below_curve1, xvals_above_curve1,
#      yvals_above_curve1, stnbelow1, stnabove1) = fit_curve_get_vals_below_curve(
#         x=x1_gd_corr, y=y1_gd_corr, func=func, stns=s1_gd_corr, shift_per=0.02)
# (y_fitted_shifted2, xvals_below_curve2,
#  yvals_below_curve2, xvals_above_curve2,
#  yvals_above_curve2, stnbelow2, stnabove2) = fit_curve_get_vals_below_curve(
#     x=x2_gd_corr, y=y2_gd_corr, func=func, stns=s2_gd_corr, shift_per=0.25)

# (y_fitted_shifted3, xvals_below_curve3,
#  yvals_below_curve3, xvals_above_curve3,
#  yvals_above_curve3, stnbelow3, stnabove3) = fit_curve_get_vals_below_curve(
#     x=x3_gd_corr, y=y3_gd_corr, func=func, stns=s3_gd_corr, shift_per=0.25)
#
# (y_fitted_shifted4, xvals_below_curve4,
#  yvals_below_curve4, xvals_above_curve4,
#  yvals_above_curve4, stnbelow4, stnabove4) = fit_curve_get_vals_below_curve(
#     x=x4_gd_corr, y=y4_gd_corr, func=func, stns=s4_gd_corr, shift_per=0.25)
#
# (y_fitted_shifted5, xvals_below_curve5,
#  yvals_below_curve5, xvals_above_curve5,
#  yvals_above_curve5, stnbelow5, stnabove5) = fit_curve_get_vals_below_curve(
#     x=x5_gd_corr, y=y5_gd_corr, func=func, stns=s5_gd_corr, shift_per=0.25)
#
# (y_fitted_shifted6, xvals_below_curve6,
#  yvals_below_curve6, xvals_above_curve6,
#  yvals_above_curve6, stnbelow6, stnabove6) = fit_curve_get_vals_below_curve(
#     x=x6_gd_corr, y=y6_gd_corr, func=func, stns=s6_gd_corr, shift_per=0.25)
#
# (y_fitted_shifted7, xvals_below_curve7,
#  yvals_below_curve7, xvals_above_curve7,
#  yvals_above_curve7, stnbelow7, stnabove7) = fit_curve_get_vals_below_curve(
#     x=x7_gd_corr, y=y7_gd_corr, func=func, stns=s7_gd_corr, shift_per=0.25)
# =============================================================================


# stns_keep_all = np.intersect1d(stnabove0, stnabove1)

# s0_keep_initial = in_df0.index.intersection(stns_keep_all)
# s1_keep_initial = in_df1.index.intersection(stns_keep_all)
# s2_keep_initial = in_df2.index.intersection(stns_keep_all)
# s3_keep_initial = in_df3.index.intersection(stns_keep_all)
# s4_keep_initial = in_df4.index.intersection(stns_keep_all)
# s5_keep_initial = in_df5.index.intersection(stns_keep_all)
# s6_keep_initial = in_df6.index.intersection(stns_keep_all)
# s7_keep_initial = in_df7.index.intersection(stns_keep_all)

# # %%
# x0_abv = in_df0.loc[s0_keep_initial, 'Distance to neighbor'].dropna().values
# y0_abv = in_df0.loc[s0_keep_initial,
#                     'Bool_Spearman_Correlation'].dropna().values

# x1_abv = in_df1.loc[s1_keep_initial, 'Distance to neighbor'].dropna().values
# y1_abv = in_df1.loc[s1_keep_initial,
#                     'Bool_Spearman_Correlation'].dropna().values

# x2_abv = in_df2.loc[s2_keep_initial, 'Distance to neighbor'].dropna().values
# y2_abv = in_df2.loc[s2_keep_initial,
#                     'Bool_Spearman_Correlation'].dropna().values

# x3_abv = in_df3.loc[s3_keep_initial, 'Distance to neighbor'].dropna().values
# y3_abv = in_df3.loc[s3_keep_initial,
#                     'Bool_Spearman_Correlation'].dropna().values
#
# x4_abv = in_df4.loc[s4_keep_initial, 'Distance to neighbor'].dropna().values
# y4_abv = in_df4.loc[s4_keep_initial,
#                     'Bool_Spearman_Correlation'].dropna().values
#
# x5_abv = in_df5.loc[s5_keep_initial, 'Distance to neighbor'].dropna().values
# y5_abv = in_df5.loc[s5_keep_initial,
#                     'Bool_Spearman_Correlation'].dropna().values
#
# x6_abv = in_df6.loc[s6_keep_initial, 'Distance to neighbor'].dropna().values
# y6_abv = in_df6.loc[s6_keep_initial,
#                     'Bool_Spearman_Correlation'].dropna().values
#
# x7_abv = in_df7.loc[s7_keep_initial, 'Distance to neighbor'].dropna().values
# y7_abv = in_df7.loc[s7_keep_initial,
#                     'Bool_Spearman_Correlation'].dropna().values


#==============================================================================
# Upper Limit
#==============================================================================
# remove_upper_limit = False
# if remove_upper_limit:
#     shit_curve_upp0 = -0.16
#     shit_curve_upp1 = -0.16
#     shit_curve_upp2 = -0.16
#
#     (y_fitted_shifted0_, xvals_below_curve0_,
#      yvals_below_curve0_, x0_abv_22,
#      y0_abv_22, stn0_below_22, stnabove0_) = fit_curve_get_vals_below_curve(
#         x=x0_abv, y=y0_abv, func=func, stns=s0_keep_initial, shift_per=shit_curve_upp0)
#     print(stn0_below_22.shape, stnabove0_.shape)

#     (y_fitted_shifted1_, xvals_below_curve1_,
#      yvals_below_curve1_, x1_abv_22,
#      y1_abv_22, stn1_below_22, stnabove1_) = fit_curve_get_vals_below_curve(
#         x=x1_abv, y=y1_abv, func=func, stns=s1_keep_initial, shift_per=shit_curve_upp1)
#     print(stn1_below_22.shape, stnabove1_.shape)

#     (y_fitted_shifted2_, xvals_below_curve2_,
#      yvals_below_curve2_, x2_abv_22,
#      y2_abv_22, stn2_below_22, stnabove2_) = fit_curve_get_vals_below_curve(
#         x=x2_abv, y=y2_abv, func=func, stns=s2_keep_initial, shift_per=shit_curve_upp1)
#
#     print(stn2_below_22.shape, stnabove2_.shape)

#     (y_fitted_shifted3_, xvals_below_curve3_,
#      yvals_below_curve3_, x3_abv_22,
#      y3_abv_22, stn3_below_22, stnabove3_) = fit_curve_get_vals_below_curve(
#         x=x3_abv, y=y3_abv, func=func, stns=s3_keep_initial, shift_per=shit_curve_upp2)
#     print(stn3_below_22.shape, stnabove3_.shape)
#
#     (y_fitted_shifted4_, xvals_below_curve4_,
#      yvals_below_curve4_, x4_abv_22,
#      y4_abv_22, stn4_below_22, stnabove4_) = fit_curve_get_vals_below_curve(
#         x=x4_abv, y=y4_abv, func=func, stns=s4_keep_initial, shift_per=shit_curve_upp2)
#     print(stn4_below_22.shape, stnabove4_.shape)
#
#     (y_fitted_shifted5_, xvals_below_curve5_,
#      yvals_below_curve5_, x5_abv_22,
#      y5_abv_22, stn5_below_22, stnabove5_) = fit_curve_get_vals_below_curve(
#         x=x5_abv, y=y5_abv, func=func, stns=s5_keep_initial, shift_per=shit_curve_upp2)
#     print(stn5_below_22.shape, stnabove5_.shape)
#
#     (y_fitted_shifted6_, xvals_below_curve6_,
#      yvals_below_curve6_, x6_abv_22,
#      y6_abv_22, stn6_below_22, stnabove6_) = fit_curve_get_vals_below_curve(
#         x=x6_abv, y=y6_abv, func=func, stns=s6_keep_initial, shift_per=shit_curve_upp2)
#     print(stn6_below_22.shape, stnabove6_.shape)
#
#     (y_fitted_shifted7_, xvals_below_curve7_,
#      yvals_below_curve7_, x7_abv_22,
#      y7_abv_22, stn7_below_22, stnabove7_) = fit_curve_get_vals_below_curve(
#         x=x7_abv, y=y7_abv, func=func, stns=s7_keep_initial, shift_per=shit_curve_upp2)
#     print(stn7_below_22.shape, stnabove7_.shape)


# #==============================================================================
# # Lower Limit
# #==============================================================================

# if percent == '95_0':
#     lim1 = 0.5
#     lim2 = 0.5
#
# if percent == '97_0':
#     lim1 = 0.45
#     lim2 = 0.42
#
# if percent == '98_0':
#     lim1 = 0.42
#     lim2 = 0.4
# #
# if percent == '99_0':
#     lim1 = 0.3
#     lim2 = 0.3
#
# x0_abv_2 = x0_abv
# y0_abv_2 = y0_abv
#
# y1_abv_2 = y1_abv[y1_abv >= lim1]
# x1_abv_2 = x1_abv[np.where(y1_abv_2 >= lim1)]
# stn1_abv_2 = s1_keep_initial[np.where(y1_abv_2 >= lim1)]
#
# y2_abv_2 = y2_abv[y2_abv >= lim1]
# x2_abv_2 = x2_abv[np.where(y2_abv_2 >= lim1)]
# stn2_abv_2 = s2_keep_initial[np.where(y2_abv_2 >= lim1)]

# y3_abv_2 = yvals_below_curve3_[yvals_below_curve3_ >= lim1]
# x3_abv_2 = xvals_below_curve3_[np.where(y3_abv_2 >= lim1)]
# stn3_abv_2 = stn3_below_22[np.where(y3_abv_2 >= lim1)]
#
# y4_abv_2 = yvals_below_curve4_[yvals_below_curve4_ >= lim2]
# x4_abv_2 = xvals_below_curve4_[np.where(y4_abv_2 >= lim2)]
# stn4_abv_2 = stn4_below_22[np.where(y4_abv_2 >= lim2)]
#
# y5_abv_2 = yvals_below_curve5_[yvals_below_curve5_ >= lim2]
# x5_abv_2 = xvals_below_curve5_[np.where(y5_abv_2 >= lim2)]
# stn5_abv_2 = stn5_below_22[np.where(y5_abv_2 >= lim2)]
#
# y6_abv_2 = yvals_below_curve6_[yvals_below_curve6_ >= lim3]
# x6_abv_2 = xvals_below_curve6_[np.where(y6_abv_2 >= lim3)]
# stn6_abv_2 = stn6_below_22[np.where(y6_abv_2 >= lim3)]
#
# y7_abv_2 = yvals_below_curve7_[yvals_below_curve7_ >= lim3]
# x7_abv_2 = xvals_below_curve7_[np.where(y7_abv_2 >= lim3)]
# stn7_abv_2 = stn7_below_22[np.where(y7_abv_2 >= lim3)]

#==============================================================================
#
#==============================================================================
# stns_keep_all_final = np.intersect1d(stnabove0, s1_keep_initial)
#                                      stn2_abv_2)
#         ,
#         stn2_abv_2), stn3_abv_2), stn4_abv_2),
#     stn5_abv_2), stn6_abv_2), stn7_abv_2)

# s01 = np.intersect1d(stn0_below_22, stn1_below_22)
# print(s01.shape)
# s012 =  np.intersect1d(s01, stn2_below_22)
# print(stns_keep_all_final.shape)
#
# s0123 =  np.intersect1d(s012, stn3_below_22)

#==============================================================================
# GET FINAL STATIONS
#==============================================================================


x0_abv_2 = in_df0.loc[stnabove0,
                      'Distance to neighbor'].dropna().values
y0_abv_2 = in_df0.loc[stnabove0,
                      'Bool_Spearman_Correlation'].dropna().values

# x1_abv_2 = in_df1.loc[stnabove0,
#                       'Distance to neighbor'].dropna().values
# y1_abv_2 = in_df1.loc[stnabove0,
#                       'Bool_Spearman_Correlation'].dropna().values

# x2_abv_2 = in_df2.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y2_abv_2 = in_df2.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x3_abv_2 = in_df3.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y3_abv_2 = in_df3.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x4_abv_2 = in_df4.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y4_abv_2 = in_df4.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x5_abv_2 = in_df5.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y5_abv_2 = in_df5.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x6_abv_2 = in_df6.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y6_abv_2 = in_df6.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x7_abv_2 = in_df7.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y7_abv_2 = in_df7.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values

#==============================================================================
#
#==============================================================================
# print('refilter')


# y0_abv_2 = y0_abv_2[y0_abv_2 >= y0_dwd.min()]
# x0_abv_2 = x0_abv[np.where(y0_abv_2 >= y0_dwd.min())]
#
# y1_abv_2 = y1_abv_2[y1_abv_2 >= y1_dwd.min()]
# x1_abv_2 = x1_abv[np.where(y1_abv_2 >= y1_dwd.min())]


# stns_keep_all_final0 = stns_keep_all_final[np.where(y0_abv_2 >= y0_dwd.min())]
# stns_keep_all_final1 = stns_keep_all_final[np.where(y1_abv_2 >= y1_dwd.min())]
# stns_keep_all_final2 = stns_keep_all_final[np.where(y2_abv_2 >= lim1)]
# stns_keep_all_final3 = stns_keep_all_final[np.where(y3_abv_2 >= lim1)]
# stns_keep_all_final4 = stns_keep_all_final[np.where(y4_abv_2 >= lim2)]
# stns_keep_all_final5 = stns_keep_all_final[np.where(y5_abv_2 >= lim2)]
# stns_keep_all_final6 = stns_keep_all_final[np.where(y6_abv_2 >= lim3)]
# stns_keep_all_final7 = stns_keep_all_final[np.where(y7_abv_2 >= lim3)]

# stns_keep_all_final = np.intersect1d(
#     stns_keep_all_final0, stns_keep_all_final1)
# #     stns_keep_all_final2)
# #         , stns_keep_all_final3), stns_keep_all_final4),
# #     stns_keep_all_final5), stns_keep_all_final6), stns_keep_all_final7)
#
# print(stns_keep_all_final.shape)


#==============================================================================
# reget all coordinates and corelations
#==============================================================================
# print(stns_keep_all_final.shape)
# x0_abv_2 = in_df0.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y0_abv_2 = in_df0.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x1_abv_2 = in_df1.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y1_abv_2 = in_df1.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values

# x1 = x1_abv_2[np.where(0 < x1_abv_2) and np.where(x1_abv_2 <= 20000)]
# s1 = s1[np.where(0 < x1_abv_2) and np.where(x1_abv_2 <= 20000)]
# y1 = y1[np.where(0 < x1_abv_2) and np.where(x1_abv_2 <= 20000)]
# in_df1 = in_df1.loc[in_df1.index.intersection(s1), :]

# stns_keep_all_final = np.intersect1d(
#     stns_keep_all_final, s1)

# idx_remove0 = []
# for ix0, yv0 in enumerate(y0_abv_2):
#
#     if yv0 <= y0_dwd.min():
#         # print(ix0)
#         idx_remove0.append(ix0)
# #         try:
#
# x0_abv_2_new = np.delete(x0_abv_2, idx_remove0)
# y0_abv_2_new = np.delete(y0_abv_2, idx_remove0)
# stns_keep_all_final_new = np.delete(stns_keep_all_final, idx_remove0)
# # #
# # #         except Exception as msg:
# x1_abv_2 = in_df1.loc[stns_keep_all_final_new,
#                       'Distance to neighbor'].dropna().values
# y1_abv_2 = in_df1.loc[stns_keep_all_final_new,
#                       'Bool_Spearman_Correlation'].dropna().values
# # #
# idx_remove1 = []
# for ix, yv1 in enumerate(y1_abv_2):
#
#     if yv1 <= y1_dwd.min():
#         # print(ix)
#         idx_remove1.append(ix)
# # #         try:
# x1_abv_2_new = np.delete(x1_abv_2, idx_remove1)
# y1_abv_2_new = np.delete(y1_abv_2, idx_remove1)
# stns_keep_all_final_new = np.delete(stns_keep_all_final_new, idx_remove1)
# #         except Exception as msg:
#             print(msg)
#             continue
# print(stns_keep_all_final_new.shape)
# #
# x0_abv_2 = in_df0.loc[stns_keep_all_final_new,
#                       'Distance to neighbor'].dropna().values
# y0_abv_2 = in_df0.loc[stns_keep_all_final_new,
#                       'Bool_Spearman_Correlation'].dropna().values
# print(stns_keep_all_final_new.shape[0])
# # stns_keep_all_final = stn0_below_22

stns_keep_al_sr = pd.DataFrame(data=stnabove0,
                               columns=['Stations'])

stns_keep_al_sr.to_csv(
    (save_dir /
        (r'keep_stns_%s_per_%s_s0_shit_%dperc.csv'
         % (percent, time_freq, shift_by_percent))),
    sep=';')

# x0_abv_2_new = in_df0.loc[
#     in_df0.index.intersection(stns_keep_all_final_new),
#     'Distance to neighbor'].dropna().values
# y0_abv_2_new = in_df0.loc[
#     in_df0.index.intersection(stns_keep_all_final_new),
#     'Bool_Spearman_Correlation'].dropna().values
# #
# x1_abv_2_new = in_df1.loc[
#     in_df1.index.intersection(stns_keep_all_final_new),
#     'Distance to neighbor'].dropna().values
# y1_abv_2_new = in_df1.loc[in_df1.index.intersection(stns_keep_all_final_new),
#                           'Bool_Spearman_Correlation'].dropna().values
#
# x2_abv_2 = in_df2.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y2_abv_2 = in_df2.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values

# x3_abv_2 = in_df3.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y3_abv_2 = in_df3.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x4_abv_2 = in_df4.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y4_abv_2 = in_df4.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x5_abv_2 = in_df5.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y5_abv_2 = in_df5.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x6_abv_2 = in_df6.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y6_abv_2 = in_df6.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x7_abv_2 = in_df7.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y7_abv_2 = in_df7.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values

#==============================================================================
#
#==============================================================================

# plt.figure(figsize=(36, 18), dpi=300)
#
# marker_size_abv_curve = 100
# marker_size_below_curve = 95
# marker_size_curve = 55

max_x = x0.max()
#
plt.ioff()

fig, axs = plt.subplots(1, 2, sharex=True, sharey=True,
                        figsize=(24, 10), dpi=300)

# plt.figure(figsize=(12, 8), dpi=300)

axs[0].scatter(x0, y0, c='r', alpha=0.75, marker='x', s=34)
# axs[0].scatter(x1_abv_2_new, y1_abv_2_new, c='b', alpha=0.75, marker='x', s=34)
# axs[0].scatter(x2_abv_2, y2_abv_2, c='g', alpha=0.75, marker='x', s=34)

# plt.scatter(x1_abv_2_new, y1_abv_2_new, c='b', alpha=0.5,
#             marker='.', label=('Second Neighbor Stn nbr %d / %d'
#                                % (y1_abv_2_new.shape[0], y1.shape[0])),
#             s=marker_size_abv_curve)
# plt.scatter(x2_abv_2, y2_abv_2, c='g', alpha=0.5,
#             marker='d', label=('Third Neighbor Stn nbr %d / %d'
#                                % (y2_abv_2.shape[0], y2.shape[0])),
#             s=marker_size_abv_curve)


#label='First Neighbor %d Pairs' % y0.shape[0],
axs[1].scatter(x0, y_fitted_shifted0, c='b', alpha=0.25,
               marker='_', s=20)

axs[1].scatter(x0_dwd, y0_dwd, c='k', alpha=0.5,
               marker='o', s=34)
axs[1].scatter(xvals_above_curve0, yvals_above_curve0,
               c='r', alpha=0.75, marker='x', s=34)
# axs[1].scatter(x1_after, y1_after, c='b', alpha=0.75, marker='x', s=34)
# axs[1].scatter(x2_after, y2_after, c='g', alpha=0.75, marker='x', s=34)


# axs[1].scatter(x0_dwd0, y0_dwd0, c='k', alpha=0.5, marker='o', s=34)
# axs[1].scatter(x1_dwd0, y1_dwd0, c='k', alpha=0.5, marker='o', s=34)
# axs[1].scatter(x2_dwd0, y2_dwd0, c='k', alpha=0.5, marker='o', s=34)


axs[0].legend(title='%d pairs' % y0.shape[0], loc='upper right',  # a) left
              frameon=False, fontsize=26)._legend_box.align = 'right'
axs[1].legend(title='%d pairs' % yvals_above_curve0.shape[0], loc='upper right',
              frameon=False, fontsize=26)._legend_box.align = 'right'


axs[0].set_xlim([0, max_x + 300])
axs[1].set_xlim([0, max_x + 300])

axs[0].set_xticks(np.arange(0, max_x + 300, 5000))
axs[1].set_xticks(np.arange(0, max_x + 300, 5000))

axs[0].set_ylim([-0.1, 1.1])
axs[0].set_xlabel('Distance [m]', labelpad=14)
axs[1].set_xlabel('Distance [m]', labelpad=14)
axs[0].set_ylabel('Indicator Correlation', labelpad=16)
# plt.legend(loc=0)
axs[0].grid(alpha=.25)
axs[1].grid(alpha=.25)
plt.tight_layout()

# plt.title('Keeping %s %d stations: Indicator correlation with distance'
#           ' for upper %s percent for %s data values'
#           % (data_source0,  x0_abv_2.shape[0],
#               percent, time_freq))

plt.savefig(save_dir /  # filtered_2_neighbors_data
            (r'_%s_%s_%s_fit_fct_%s_shit_by_%dperc.png'
             % (data_source0, data_source, percent, time_freq, shift_by_percent)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
