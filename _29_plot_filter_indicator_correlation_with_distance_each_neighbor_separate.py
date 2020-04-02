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
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from pathlib import Path
from math import exp
from _00_additional_functions import (list_all_full_path,
                                      fit_exp_fct_get_vals_below_abv,
                                      fit_curve_get_vals_below_curve,
                                      gen_path_df_file,
                                      read_filter_df_corr_return_stns_x_y_vals,
                                      remove_all_low_corr_short_dist)

plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'axes.labelsize': 26})

#==============================================================================
#
#==============================================================================
# main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

# for BW
# data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_DWD_ppt_correlation_'
# data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'

# for RH
# data_dir_Netamto_dfs = Path(
#    r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\indicator_correlation')

data_dir_Netamto_dfs = (r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz/indicator_correlation')


data_dir_DWD_dfs = data_dir_Netamto_dfs
netatmo_path_acc_b4 = r'1pearson_year_allyears_df_comparing_correlations_max_sep_dist_1000000_'

netatmo_path_acc = r'pearson_year_allyears_df_comparing_correlations_max_sep_dist_1000000_'

dwd_path_Acc = r'pearson_year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '99_0' 
time_freq = '60min'

data_source0 = 'Netatmo'  # reference station 'Netatmo'

data_source = 'dwd'  # compare to station 'netatmo'

remove_upper_limit = False

plot_dwd_on_top = True

shift_factor = 20
shift_by_percent = 10
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
    in_df0_dwd = pd.read_csv(df0_dwd, index_col=0, sep=';').dropna(how='all')
    in_df1_dwd = pd.read_csv(df1_dwd, index_col=0, sep=';').dropna(how='all')
    in_df2_dwd = pd.read_csv(df2_dwd, index_col=0, sep=';').dropna(how='all')
    x0_dwd = in_df0_dwd.loc[:, 'Distance to neighbor'].values.ravel()
    x1_dwd = in_df1_dwd.loc[:, 'Distance to neighbor'].values.ravel()
    x2_dwd = in_df2_dwd.loc[:, 'Distance to neighbor'].values.ravel()
    
    # x3_dwd0 = in_df3_dwd.loc[:, 'Distance to neighbor'].values.ravel()
    # x4_dwd0 = in_df4_dwd.loc[:, 'Distance to neighbor'].values.ravel()
    # x5_dwd0 = in_df5_dwd.loc[:, 'Distance to neighbor'].values.ravel()
    
    y0_dwd = in_df0_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
    y1_dwd = in_df1_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
    y2_dwd = in_df2_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
#==============================================================================

# %% apply filter for every neighbor seperatly

def exp_func(x, b): 
    return  np.exp(-x*b)

x = np.linspace(0, max(x0), x0.size)
x_scaled = x / max(x0)

"""
popt, pcov = curve_fit(exp_func, x_scaled, y0)

plt.ioff()
plt.figure()
plt.scatter(x0, y0, c='r', marker='x', label="Original Data")

plt.plot(x, exp_func(x_scaled, *popt), 'y-', label="Fitted Curve")
plt.plot(x, exp_func(x_scaled, *popt*1.5), c='b', label="Fitted Curve")

plt.plot(x, exp_func(x_scaled, *popt*2), 'c-', label="Fitted Curve")
plt.plot(x, exp_func(x_scaled, *popt*3), 'm-', label="Fitted Curve")

plt.plot(x, exp_func(x_scaled, *popt*4), 'orange', label="Fitted Curve")
plt.plot(x, exp_func(x_scaled, *popt*6), 'g-', label="Fitted Curve")
plt.plot(x, exp_func(x_scaled, *popt*7), 'lightblue', label="Fitted Curve")

plt.plot(x, exp_func(x_scaled, *popt*8), 'grey', label="Fitted Curve")
plt.plot(x, exp_func(x_scaled, *popt*9), 'k', label="Fitted Curve")
plt.plot(x, exp_func(x_scaled, *popt*12), 'k', label="Fitted Curve")
plt.grid(alpha=.5)
plt.show()
#plt.legend()
"""

def fit_exp_fct_get_vals_below_abv(x, y, exp_func, stns, shift_factor=1):
    """ fit exponential function to data and shifted 10% downwards
    return all data above and below function with 
    corresponding distance and correlation values and station ids
    """

    # bounds=[[a1,b1],[a2,b2]]
    #x = np.insert(x, 0,0)
    #y = np.insert(y, 1,0)
    
    # x_arr = np.linspace(0, max(x), x.size)
    x_scaled = x / max(x)
    
    popt, _ = curve_fit(exp_func, x_scaled, y)
    
    
    print('fitted parameters are %.2f' % popt)  # , popt[3])
    #y_fitted = exp_func(x_scaled, *popt)
    
    #exp_fitted = exp_func(x_scaled, *popt)
    
    shifted_param = popt * shift_factor
    y_fitted_shifted = exp_func(x_scaled, *shifted_param)

    xvals_below_curve = x[np.where(y <= y_fitted_shifted)]
    yvals_below_curve = y[np.where(y <= y_fitted_shifted)]
    
    stns_below_curve = stns[np.where(y <= y_fitted_shifted)]

    xvals_above_curve = x[np.where(y > y_fitted_shifted)]
    yvals_above_curve = y[np.where(y > y_fitted_shifted)]
        
    stns_above_curve = stns[np.where(y > y_fitted_shifted)]
    
    #test_ = in_df0.loc[stns_above_curve, 'Distance to neighbor'].values
    #test_2 = in_df0.loc[stns_above_curve, 'Bool_Spearman_Correlation'].values
    
    #plt.ioff()
    #plt.scatter(xvals_below_curve, yvals_below_curve, c='r')
    #plt.scatter(test_, test_2, c='g')
    #plt.scatter(x0, y_fitted_shifted, c='k')
    #plt.scatter(xvals_above_curve, yvals_above_curve, c='b')
    #plt.show()
        # plt.scatter(x_arr, exp_fitted)
    # 

    
    return (y_fitted_shifted, xvals_below_curve, yvals_below_curve,
            xvals_above_curve, yvals_above_curve, stns_below_curve,
            stns_above_curve)
    
#shift_factor = 5
(y_fitted_shifted0_exp, xvals_below_curve0_exp,
 yvals_below_curve0_exp, xvals_above_curve0_exp,
 yvals_above_curve0_exp, stnbelow0_exp, stnabove0_exp) = fit_exp_fct_get_vals_below_abv(
    x=x0, y=y0, exp_func=exp_func, stns=s0, shift_factor=shift_factor)
 
 
#x0 = np.insert(x0, 0 ,0)
#y0 = np.insert(y0, 0, 1)
#s0 = np.insert(s0, 0, s0[0])
#shift_by_percent = 20
(y_fitted_shifted0, xvals_below_curve0,
 yvals_below_curve0, xvals_above_curve0,
 yvals_above_curve0, stnbelow0, stnabove0) = fit_curve_get_vals_below_curve(
    x=xvals_above_curve0_exp, y=yvals_above_curve0_exp,
     stns=stnabove0_exp, shift_per=shift_by_percent / 100)

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



# #==============================================================================
# # Lower Limit
# #==============================================================================

# if percent == '95_0':
#     lim1 = 0.5
#     lim2 = 0.5


# x0_abv_2 = x0_abv
# y0_abv_2 = y0_abv
#
# y1_abv_2 = y1_abv[y1_abv >= lim1]
# x1_abv_2 = x1_abv[np.where(y1_abv_2 >= lim1)]
# stn1_abv_2 = s1_keep_initial[np.where(y1_abv_2 >= lim1)]


#==============================================================================
# GET FINAL STATIONS
#==============================================================================


x0_abv_2 = in_df0.loc[stnabove0,
                      'Distance to neighbor'].dropna().values
y0_abv_2 = in_df0.loc[stnabove0,
                      'Bool_Spearman_Correlation'].dropna().values

y0_shifted = in_df0.loc[s0,
                      'Bool_Spearman_Correlation_Shifted'].dropna().values
#==============================================================================
#
#==============================================================================
# print('refilter')


# y0_abv_2 = y0_abv_2[y0_abv_2 >= y0_dwd.min()]
# x0_abv_2 = x0_abv[np.where(y0_abv_2 >= y0_dwd.min())]


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
## # stns_keep_all_final = stn0_below_22

stns_keep_al_sr = pd.DataFrame(data=stnabove0,
                               columns=['Stations'])

stns_keep_al_sr.to_csv(
    os.path.join(save_dir,
        (r'keep_stns_%s_per_%s_shift_%dperc_%dfact.csv'
         % (percent, time_freq, shift_by_percent, shift_factor))),
    sep=';')

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
# axs[0].scatter(x0, y0_shifted, c='b', alpha=0.75, marker='o', s=34)

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
axs[1].scatter(x0, y_fitted_shifted0_exp, c='darkblue', alpha=0.25,
               marker='_', s=20)

axs[1].scatter(xvals_above_curve0_exp, y_fitted_shifted0, c='darkgreen',
               alpha=0.25,
               marker='_', s=20)

axs[1].scatter(xvals_below_curve0_exp, yvals_below_curve0_exp,
               c='b', alpha=0.75, marker='.', s=34)

axs[1].scatter(xvals_below_curve0, yvals_below_curve0,
               c='g', alpha=0.75, marker='d', s=34)

axs[1].scatter(xvals_above_curve0, yvals_above_curve0,
               c='r', alpha=0.75, marker='x', s=34)

if plot_dwd_on_top:
    axs[1].scatter(x0_dwd, y0_dwd, c='k', alpha=0.5,
               marker='o', s=34)
    
    axs[1].scatter(x1_dwd, y1_dwd, c='k', alpha=0.5,
               marker='o', s=34)
    
#     axs[1].scatter(x2_dwd, y2_dwd, c='k', alpha=0.5,
#                marker='o', s=34)
    
#axs[1].scatter(x0_abv_2, y0_abv_2, c='k', alpha=0.75, marker='x', s=34)
# axs[1].scatter(x2_after, y2_after, c='g', alpha=0.75, marker='x', s=34)


# axs[1].scatter(x0_dwd0, y0_dwd0, c='k', alpha=0.5, marker='o', s=34)
# axs[1].scatter(x1_dwd0, y1_dwd0, c='k', alpha=0.5, marker='o', s=34)
# axs[1].scatter(x2_dwd0, y2_dwd0, c='k', alpha=0.5, marker='o', s=34)


axs[0].legend(title='%d pairs' % y0.shape[0], loc='upper right',  # a) left
              frameon=False, fontsize=26)._legend_box.align = 'right'
axs[1].legend(title='%d pairs' % stnabove0.shape[0], loc='upper right',
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

plt.savefig(os.path.join(save_dir,  # filtered_2_neighbors_data
            (r'4_%s_%s_%s_fit_exp_fct_%s_shift_%dperc_%dfact_4.png'
             % (data_source0, data_source,
                 percent, time_freq, shift_by_percent,
                 shift_factor))),
            papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
