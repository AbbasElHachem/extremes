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

from _00_additional_functions import fit_curve_get_vals_below_curve


def model(x, a, b, c, d):
    """ 3degree polynomial function used for fitting and as filter"""
    return a * x**3 + b * x**2 + c * x + d

# def model(x, a, b, c):
#     return (a / np.sqrt(x + b) + c)
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
dwd_stns = np.hstack([s0_dwd, s1_dwd])
#
# # Netatmo DATA
xall0 = np.hstack([x0_gd_corr, x1_gd_corr])
yall0 = np.hstack([y0_gd_corr, y1_gd_corr])
stnsall0 = np.hstack([s0_gd_corr, s1_gd_corr])
# #%%
(y_fitted_shifted6, xvals_below_curve6,
 yvals_below_curve6, xvals_above_curve6,
 yvals_above_curve6, stnbelow6, stnabove6) = fit_curve_get_vals_below_curve(
    x=xall, y=yall, func=func, stns=dwd_stns, shift_per=0.25)
#

# #xall_ = np.insert(arr=xall, obj=0, values=0)
# #yall_ = np.insert(arr=yall, obj=0, values=1)
# #%%
guesses = [-0.00001, 0.001, 0.001, 1]
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

y_netatmo_keep = yall0[yall0 >= yvals_above_curve6.min()]
x_netatmo_keep = xall0[yall0 >= yvals_above_curve6.min()]
stn_netatmo_keep = np.unique(stnsall0[yall0 >= yvals_above_curve6.min()])
print('keeping stations', stn_netatmo_keep.shape, ' from ', stnsall0.shape)
# %%
plt.ioff()
# plt.plot(xall0,
#          yall0,
#          '.', alpha=0.8, label='Test data')
plt.scatter(xall, y_fitted_shifted6, color='r',
            zorder=0.9, label='Edge', alpha=.5)
plt.scatter(xall, yall, color='b', zorder=0.9, label='DWD Data', alpha=.5)

# plt.scatter(xall0, yall0, color='g', zorder=0.9,
#             label='Netamo Data', alpha=0.5)

plt.scatter(x_netatmo_keep, y_netatmo_keep, color='m', zorder=0.9,
            label='Netamo Data', alpha=0.5)

plt.legend(loc='lower left')
plt.show()

# =============================================================================
