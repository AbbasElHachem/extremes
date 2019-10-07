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

from scipy.optimize import leastsq
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

netatmo_path_acc = r'year_allyears_df_comparing_correlations_max_sep_dist_30000_'
dwd_path_Acc = r'year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '95'
time_freq = '60min'

data_source0 = 'Netatmo'  # reference station 'Netatmo'
data_source01 = 'DWD'  # reference station 'DWD'

data_source = 'dwd'  # compare to station 'netatmo'

remove_upper_limit = False
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
if data_source01 == 'DWD':
    # for DWD stations neighbors start from 1 not 0 !
    df0_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent, 1)
    df1_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent, 2)
    df2_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                               data_source, percent, 3)
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

#==============================================================================


s0, x0, y0, in_df0 = read_filter_df_corr_return_stns_x_y_vals(df0)
s1, x1, y1, in_df1 = read_filter_df_corr_return_stns_x_y_vals(df1)

s0_dwd, x0_dwd, y0_dwd, in_df0_dwd = read_filter_df_corr_return_stns_x_y_vals(
    df0_dwd)
s1_dwd, x1_dwd, y1_dwd, in_df1_dwd = read_filter_df_corr_return_stns_x_y_vals(
    df1_dwd)
# s2, x2, y2, in_df2 = read_filter_df_corr_return_stns_x_y_vals(df2)
# s3, x3, y3, in_df3 = read_filter_df_corr_return_stns_x_y_vals(df3)
# s4, x4, y4, in_df4 = read_filter_df_corr_return_stns_x_y_vals(df4)
# s5, x5, y5, in_df5 = read_filter_df_corr_return_stns_x_y_vals(df5)
# s6, x6, y6, in_df6 = read_filter_df_corr_return_stns_x_y_vals(df6)
# s7, x7, y7, in_df7 = read_filter_df_corr_return_stns_x_y_vals(df7)
#==============================================================================

# x0, y0, s0 = remove_all_low_corr_short_dist(x0, y0, 5e3, 0.6, s0)
x0_gd_corr, y0_gd_corr, s0_gd_corr, df0_gd_corr = remove_all_low_corr_short_dist(
    x0, y0, 5e3, 0.4, s0, in_df0)

x1_gd_corr, y1_gd_corr, s1_gd_corr, df1_gd_corr = remove_all_low_corr_short_dist(
    x1, y1, 5e3, 0.4, s1, in_df1)

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


def model(x, a, b, c):
    return (a / np.sqrt(x + b) + c)
#%%


def get_flipped(y_data, y_model):
    flipped = y_model - y_data
    flipped[flipped > 0] = 0
    return flipped


def flipped_resid(pars, x, y):
    """
    For every iteration, everything above the currently proposed
    curve is going to be mirrored down, so that the next iterations
    is going to progressively shift downwards.
    """
    y_model = model(x, *pars)
    flipped = get_flipped(y, y_model)
    resid = np.square(y + flipped - y_model)
    # print pars, resid.sum() # uncomment to check the iteration parameters
    return np.nan_to_num(resid)


# # DWD DATA
xall = np.hstack([x0_dwd, x1_dwd])
yall = np.hstack([y0_dwd, y1_dwd])

#
# # Netatmo DATA
xall0 = np.hstack([x0_gd_corr, x1_gd_corr])
yall0 = np.hstack([y0_gd_corr, y1_gd_corr])
stnsall0 = np.hstack([s0_gd_corr, s1_gd_corr])
# #%%
#
# #xall_ = np.insert(arr=xall, obj=0, values=0)
# #yall_ = np.insert(arr=yall, obj=0, values=1)
# #%%
guesses = [100, 100, 0]
#
fit_pars, flag = leastsq(func=flipped_resid, x0=guesses,
                         args=(xall, yall))
#
y_fit = model(xall, *fit_pars)
#
# # shift curves 5percent downwards
y_fit = y_fit - 0.025 * y_fit.min()
#%%
guesses = fit_pars

# fit for Netatmo, not needed
y_fit0 = model(xall0, *fit_pars)

y_netatmo_keep = yall0[yall0 >= y_fit.min()]
x_netatmo_keep = xall0[yall0 >= y_fit.min()]
stn_netatmo_keep = np.unique(stnsall0[yall0 >= y_fit.min()])
print('keeping stations', stn_netatmo_keep.shape, ' from ', stnsall0.shape)
# %%
plt.ioff()
# plt.plot(xall0,
#          yall0,
#          '.', alpha=0.8, label='Test data')
plt.scatter(xall, y_fit, color='r', zorder=0.9, label='Edge', alpha=.5)
plt.scatter(xall, yall, color='b', zorder=0.9, label='DWD Data', alpha=.5)

plt.scatter(xall0, yall0, color='g', zorder=0.9,
            label='Netamo Data', alpha=0.5)

plt.scatter(x_netatmo_keep, y_netatmo_keep, color='m', zorder=0.9,
            label='Netamo Data', alpha=0.5)

plt.legend(loc='lower left')
plt.show()

# =============================================================================

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
# stns_keep_all_final = np.intersect1d(stnabove0, stn1_abv_2)
#                                      stn2_abv_2)
#         ,
#         stn2_abv_2), stn3_abv_2), stn4_abv_2),
#     stn5_abv_2), stn6_abv_2), stn7_abv_2)

# s01 = np.intersect1d(stn0_below_22, stn1_below_22)
# print(s01.shape)
# s012 =  np.intersect1d(s01, stn2_below_22)
# print(s012.shape)
#
# s0123 =  np.intersect1d(s012, stn3_below_22)

#==============================================================================
# GET FINAL STATIONS
#==============================================================================

#
# x0_abv_2 = in_df0.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y0_abv_2 = in_df0.loc[stns_keep_all_final,
#                       'Bool_Spearman_Correlation'].dropna().values
#
# x1_abv_2 = in_df1.loc[stns_keep_all_final,
#                       'Distance to neighbor'].dropna().values
# y1_abv_2 = in_df1.loc[stns_keep_all_final,
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
# lim1 = 0.5
# lim2 = 0.48
# lim3 = 0.44
#
# x0_abv_2 = x0_abv_2
# y0_abv_2 = y0_abv_2
#
# stns_keep_all_final0 = stns_keep_all_final[np.where(y0_abv_2 >= lim1)]
# stns_keep_all_final1 = stns_keep_all_final[np.where(y1_abv_2 >= lim2)]
# # stns_keep_all_final2 = stns_keep_all_final[np.where(y2_abv_2 >= lim1)]
# stns_keep_all_final3 = stns_keep_all_final[np.where(y3_abv_2 >= lim1)]
# stns_keep_all_final4 = stns_keep_all_final[np.where(y4_abv_2 >= lim2)]
# # stns_keep_all_final5 = stns_keep_all_final[np.where(y5_abv_2 >= lim2)]
# # stns_keep_all_final6 = stns_keep_all_final[np.where(y6_abv_2 >= lim3)]
# # stns_keep_all_final7 = stns_keep_all_final[np.where(y7_abv_2 >= lim3)]
#
# stns_keep_all_final = np.intersect1d(
#     stns_keep_all_final0, stns_keep_all_final1)
# #     stns_keep_all_final2)
# #         , stns_keep_all_final3), stns_keep_all_final4),
# #     stns_keep_all_final5), stns_keep_all_final6), stns_keep_all_final7)
#
# print(stns_keep_all_final.shape)
#
# # stns_keep_all_final = stn0_below_22
# stns_keep_al_sr = pd.DataFrame(data=stns_keep_all_final,
#                                columns=['Stations'])
#
# stns_keep_al_sr.to_csv(
#     (save_dir /
#         (r'keep_stns_all_neighbor_%s_per_%s_s0.csv'
#          % (percent, time_freq))),
#     sep=';')
# #==============================================================================
# # reget all coordinates and corelations
# #==============================================================================
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
#
# # x2_abv_2 = in_df2.loc[stns_keep_all_final,
# #                       'Distance to neighbor'].dropna().values
# # y2_abv_2 = in_df2.loc[stns_keep_all_final,
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
# plt.ioff()
# plt.figure(figsize=(36, 18), dpi=300)
#
# marker_size_abv_curve = 100
# marker_size_below_curve = 95
# marker_size_curve = 55
#
# plt.figure(figsize=(16, 12), dpi=300)
#
# marker_size_abv_curve = 34
# marker_size_below_curve = 28
# marker_size_curve = 15
#
# plt.scatter(x0_abv_2, y0_abv_2, c='r', alpha=0.5,
#             marker='x', label=('First Neighbor Stn nbr %d / %d'
#                                % (y0_abv_2.shape[0], y0.shape[0])),
#             s=marker_size_abv_curve)
# plt.scatter(x1_abv_2, y1_abv_2, c='b', alpha=0.5,
#             marker='.', label=('Second Neighbor Stn nbr %d / %d'
#                                % (y1_abv_2.shape[0], y1.shape[0])),
#             s=marker_size_abv_curve)
# # plt.scatter(x2_abv_2, y2_abv_2, c='g', alpha=0.5,
#             marker='d', label=('Third Neighbor Stn nbr %d / %d'
#                                % (y2_abv_2.shape[0], y2.shape[0])),
#             s=marker_size_abv_curve)
# plt.scatter(x3_abv_2, y3_abv_2, c='darkorange', alpha=0.5,
#             marker='*', label=('Fourth Neighbor Stn nbr %d / %d'
#                                % (y3_abv_2.shape[0], y3.shape[0])),
#             s=marker_size_abv_curve)
# plt.scatter(x4_abv_2, y4_abv_2, c='m', alpha=0.5,
#             marker='+', label=('Fifth Neighbor Stn nbr %d / %d'
#                                % (y4_abv_2.shape[0], y4.shape[0])),
#             s=marker_size_abv_curve)
# plt.scatter(x5_abv_2, y5_abv_2, c='k', alpha=0.5,
#             marker='1', label=('Sixth Neighbor Stn nbr %d / %d'
#                                % (y5_abv_2.shape[0], y5.shape[0])),
#             s=marker_size_abv_curve)
#
# plt.scatter(x6_abv_2, y6_abv_2, c='darkblue', alpha=0.5,
#             marker='X', label=('Seventh Neighbor Stn nbr %d / %d'
#                                % (y6_abv_2.shape[0], y6.shape[0])),
#             s=marker_size_abv_curve)
#
# plt.scatter(x7_abv_2, y7_abv_2, c='c', alpha=0.5,
#             marker='8', label=('Eigth Neighbor Stn nbr %d / %d'
#                                % (y7_abv_2.shape[0], y7.shape[0])),
#             s=marker_size_abv_curve)


# plt.scatter(xvals_below_curve0, yvals_below_curve0, c='grey', alpha=0.65,
#             marker='x', s=marker_size_below_curve,
#             label='First Neighbor Stn nbr %d' % yvals_below_curve0.shape[0])
# plt.scatter(xvals_below_curve1, yvals_below_curve1, c='grey', alpha=0.65,
#             marker='.',  s=marker_size_below_curve,
#             label='Second Neighbor Stn nbr %d' % yvals_below_curve1.shape[0])
# plt.scatter(xvals_below_curve2, yvals_below_curve2, c='grey', alpha=0.65,
#             marker='d',  s=marker_size_below_curve,
#             label='Third Neighbor Stn nbr %d' % yvals_below_curve2.shape[0])
# plt.scatter(xvals_below_curve3, yvals_below_curve3, c='grey', alpha=0.65,
#             marker='*',  s=marker_size_below_curve,
#             label='Fourth Neighbor Stn nbr %d' % yvals_below_curve3.shape[0])
# plt.scatter(xvals_below_curve4, yvals_below_curve4, c='grey', alpha=0.65,
#             marker='+', s=marker_size_below_curve,
#             label='Fifth Neighbor Stn nbr %d' % yvals_below_curve4.shape[0])
#
# plt.scatter(x0_abv, y_fitted_shifted0_, c='darkred', alpha=0.25,
#             marker='.',  s=marker_size_curve)  # label='Fitted curve 1',
# plt.scatter(x1_abv, y_fitted_shifted1_, c='darkblue', alpha=0.25,
#             marker='.', s=marker_size_curve)  # label='Fitted curve 2',
# plt.scatter(x2_abv, y_fitted_shifted2_, c='darkgreen', alpha=0.25,
#             marker='.', s=marker_size_curve)  # label='Fitted curve 3',
# plt.scatter(x3_abv, y_fitted_shifted3_, c='darkorange', alpha=0.25,
#             marker='.',  s=marker_size_curve)  # label='Fitted curve 4',
# plt.scatter(x4_abv, y_fitted_shifted4_, c='m', alpha=0.25,
#             marker='.', s=marker_size_curve)  # label='Fitted curve 5',
#
# plt.scatter(x5_abv, y_fitted_shifted5_, c='k', alpha=0.25,
#             marker='.',  s=marker_size_curve)  # label='Fitted curve 1',
# plt.scatter(x6_abv, y_fitted_shifted6_, c='darkblue', alpha=0.25,
#             marker='.', s=marker_size_curve)  # label='Fitted curve 2',
# plt.scatter(x7_abv, y_fitted_shifted7_, c='c', alpha=0.25,
#             marker='.', s=marker_size_curve)  # label='Fitted curve 3',

# plt.show()
#
# plt.xlim([0, max([x0_abv.max(), x1_gd_corr.max()]) + 1000])
# plt.xticks(np.arange(0, x1_gd_corr.max() + 1000, 5000))
# plt.ylim([-0.1, 1.1])
# plt.xlabel('Distance (m)')
# plt.ylabel('Indicator Spearman Correlation')
# plt.legend(loc=0)
# plt.grid(alpha=.25)
#
# plt.title('Keeping %s %d stations: Indicator correlation with distance'
#           ' for upper %s percent for %s data values'
#           % (data_source0, stns_keep_all_final.shape[0],
#               percent, time_freq))
# plt.savefig(save_dir /
#             (r'_%s_%s_%s_percent_indic_corr_freq_%s_filtered_2_neighbors_data_.png'
#              % (data_source0, data_source, percent, time_freq)),
#             frameon=True, papertype='a4',
#             bbox_inches='tight', pad_inches=.2)
# plt.close()
# #
# #def model(x, a, b, c, d):
# #    """ 3degree polynomial function used for fitting and as filter"""
# #    return a * x**3 + b * x**2 + c * x + d
#

#
#
#
# #%%
# #def get_x_y_vals_df(df, xcol='Distance to neighbor')
# # DWD data
# xvals0 = in_df0.loc[:, 'Distance to neighbor'].values.ravel()
# yvals0 = in_df0.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# xvals1 = in_df1.loc[:, 'Distance to neighbor'].values.ravel()
# yvals1 = in_df1.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# xvals2 = in_df2.loc[:, 'Distance to neighbor'].values.ravel()
# yvals2 = in_df2.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# xvals3 = in_df3.loc[:, 'Distance to neighbor'].values.ravel()
# yvals3 = in_df3.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# xvals4 = in_df4.loc[:, 'Distance to neighbor'].values.ravel()
# yvals4 = in_df4.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# #%%
# # Netatmo DATA
#
# #def get_x_y_vals_df(df, xcol='Distance to neighbor')
# xvals00 = in_df00.loc[:, 'Distance to neighbor'].values.ravel()
# yvals00 = in_df00.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# xvals01 = in_df01.loc[:, 'Distance to neighbor'].values.ravel()
# yvals01 = in_df01.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# xvals02 = in_df02.loc[:, 'Distance to neighbor'].values.ravel()
# yvals02 = in_df02.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# xvals03 = in_df03.loc[:, 'Distance to neighbor'].values.ravel()
# yvals03 = in_df03.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# xvals04 = in_df04.loc[:, 'Distance to neighbor'].values.ravel()
# yvals04 = in_df04.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#
# #%%
# # DWD DATA
# xall = np.hstack([xvals0, xvals1, xvals2, xvals3, xvals4])
# yall = np.hstack([yvals0, yvals1, yvals2, yvals3, yvals4])
#
# # Netatmo DATA
# xall0 = np.hstack([xvals00, xvals01, xvals02, xvals03, xvals04])
# yall0 = np.hstack([yvals00, yvals01, yvals02, yvals03, yvals04])
# #%%
#
# #xall_ = np.insert(arr=xall, obj=0, values=0)
# #yall_ = np.insert(arr=yall, obj=0, values=1)
# #%%
# guesses =[100, 100, 0]
#
# fit_pars, flag = leastsq(func = flipped_resid, x0 = guesses,
#                          args = (xall, yall))
#
# y_fit = model(xall, *fit_pars)
#
# # shift curves 5percent downwards
# y_fit = y_fit - 0.025*y_fit.min()
# #%%
# guesses = fit_pars
#
# fit_pars0, flag0 = leastsq(func = flipped_resid, x0 = guesses,
#                          args = (xall0, yall0))
#
# y_fit0 = model(xall0, *fit_pars)
#
# # shift curves 5percent downwards
# y_fit0 = y_fit0 - 0.025*y_fit0.min()
# # %%
# plt.ioff()
# plt.plot(xall0,
#          yall0,
#          '.', alpha=0.8, label = 'Test data')
# plt.scatter(xall, y_fit, color='r', zorder = 0.9, label = 'Edge')
# plt.scatter(xall0, y_fit0, color='g', zorder = 0.9, label = 'Guess')
# plt.legend(loc = 'lower left')
# plt.show()
# #%%

# ybelow = yall0[yall0 < y_fit0]
