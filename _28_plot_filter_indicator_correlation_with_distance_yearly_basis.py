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
import math
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
# data_dir_Netamto_dfs = (r'/run/media/abbas/EL Hachem 2019/home_office/indicator_correlation_BW_new')

# for RH
# data_dir_Netamto_dfs = Path(
#    r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\indicator_correlation')

data_dir_Netamto_dfs = (r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz/indicator_correlation')

data_dir_DWD_dfs = data_dir_Netamto_dfs
netatmo_path_acc_b4 = r'0pearson_yearly_df_comparing_correlations_max_sep_dist_1000000_'

netatmo_path_acc = r'pearson_year_allyears_df_comparing_correlations_max_sep_dist_1000000_'

dwd_path_Acc = r'pearson_year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '99_0' 
time_freq = '60min'

data_source0 = 'Netatmo'  # reference station 'Netatmo'

data_source = 'dwd'  # compare to station 'netatmo'

remove_upper_limit = False

plot_dwd_on_top = True

shift_factor = 10
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
# df1_b4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_b4, time_freq,
#                           data_source, percent, 1)
# df2_b4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_b4, time_freq,
#                           data_source, percent, 2)
# 
# df0_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 1)
# df1_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 2)
# df2_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 3)

save_dir = data_dir_Netamto_dfs
#==============================================================================

# in_df0_dwd = pd.read_csv(df0_dwd, index_col=0, sep=';').dropna(how='all')
# in_df1_dwd = pd.read_csv(df1_dwd, index_col=0, sep=';').dropna(how='all')
# in_df2_dwd = pd.read_csv(df2_dwd, index_col=0, sep=';').dropna(how='all')

# s0_dwd, x0_dwd, y0_dwd, in_df0_dwd = read_filter_df_corr_return_stns_x_y_vals(
#     df0_dwd)
# s1_dwd, x1_dwd, y1_dwd, in_df1_dwd = read_filter_df_corr_return_stns_x_y_vals(
#     df1_dwd)
# #==============================================================================
# 
#==============================================================================
# read dataframe of Netatmo
in_df0 = pd.read_csv(df0_b4, sep=';', index_col=0)

# get all data
dict_year_dist_corr = {year:{'stns': [],
                             'dist': [],
                             'corr': []}
                        for year in in_df0.index}
for year in in_df0.index:
    for stn in in_df0.loc[year, :].index:
        try:
            if math.isnan(in_df0.loc[year, stn]):
                pass
        except TypeError:
            dist_corr = in_df0.loc[year, stn].split(",")
            dist = float(dist_corr[0].replace('[', ''))
            corr = float(dist_corr[1].replace(']', ''))
            dict_year_dist_corr[year]['stns'].append(stn)
            dict_year_dist_corr[year]['dist'].append(dist)
            dict_year_dist_corr[year]['corr'].append(corr)

#==============================================================================

# %% apply filter for every neighbor seperatly


def exp_func(x, b): 
    return  np.exp(-x * b)

def fit_exp_fct_get_vals_below_abv(x, y, exp_func, stns, shift_factor=1):
    """ fit exponential function to data and shifted 10% downwards
    return all data above and below function with 
    corresponding distance and correlation values and station ids
    """

    x_scaled = x / max(x)
    
    popt, _ = curve_fit(exp_func, x_scaled, y)
    
    print('fitted parameters are %.2f' % popt)  # , popt[3])
    
    shifted_param = popt * shift_factor
    y_fitted_shifted = exp_func(x_scaled, *shifted_param)

    xvals_below_curve = x[np.where(y <= y_fitted_shifted)]
    yvals_below_curve = y[np.where(y <= y_fitted_shifted)]
    
    stns_below_curve = stns[np.where(y <= y_fitted_shifted)]

    xvals_above_curve = x[np.where(y > y_fitted_shifted)]
    yvals_above_curve = y[np.where(y > y_fitted_shifted)]
        
    stns_above_curve = stns[np.where(y > y_fitted_shifted)]
    
    return (y_fitted_shifted, xvals_below_curve, yvals_below_curve,
            xvals_above_curve, yvals_above_curve, stns_below_curve,
            stns_above_curve)
    
# # shift_factor = 5
# (y_fitted_shifted0_exp, xvals_below_curve0_exp,
#  yvals_below_curve0_exp, xvals_above_curve0_exp,
#  yvals_above_curve0_exp, stnbelow0_exp, stnabove0_exp) = fit_exp_fct_get_vals_below_abv(
#     x=x0, y=y0, exp_func=exp_func, stns=s0, shift_factor=shift_factor)

# # shift_by_percent = 20
# (y_fitted_shifted0, xvals_below_curve0,
#  yvals_below_curve0, xvals_above_curve0,
#  yvals_above_curve0, stnbelow0, stnabove0) = fit_curve_get_vals_below_curve(
#     x=xvals_above_curve0_exp, y=yvals_above_curve0_exp,
#      stns=stnabove0_exp, shift_per=shift_by_percent / 100)

#==============================================================================
# GET FINAL STATIONS
#==============================================================================
# 
# x0_abv_2 = in_df0.loc[stnabove0,
#                       'Distance to neighbor'].dropna().values
# y0_abv_2 = in_df0.loc[stnabove0,
#                       'Bool_Spearman_Correlation'].dropna().values

#==============================================================================
#
#==============================================================================
# # # stns_keep_all_final = stn0_below_22

# stns_keep_al_sr = pd.DataFrame(data=stnabove0,
#                                columns=['Stations'])
# 
# stns_keep_al_sr.to_csv(
#     os.path.join(save_dir,
#         (r'keep_stns_%s_per_%s_shift_%dperc_%dfact.csv'
#          % (percent, time_freq, shift_by_percent, shift_factor))),
#     sep=';')

#==============================================================================
#
#==============================================================================


plt.ioff()

fig, axs = plt.subplots(1, 1, sharex=True, sharey=True,
                        figsize=(12, 8), dpi=300)
colors = ['r', 'b', 'g', 'orange', 'm']
markers = ['x', 'o', ',', '.', '2']
for ix, year in enumerate(dict_year_dist_corr.keys()):
    dist = dict_year_dist_corr[year]['dist']
    corr = dict_year_dist_corr[year]['corr']
    
    axs.scatter(dist, corr, c=colors[ix], alpha=0.75,
                   marker=markers[ix], s=40)  # label=str(year))

plt.legend(loc=0)
plt.grid(alpha=.25)

plt.tight_layout()

plt.xlim([0, 2.01e4])
plt.ylim([-0.05, 1.05])
# plt.title('Keeping %s %d stations: Indicator correlation with distance'
#           ' for upper %s percent for %s data values'
#           % (data_source0,  x0_abv_2.shape[0],
#               percent, time_freq))

# plt.savefig(os.path.join(save_dir,  # filtered_2_neighbors_data
#             (r'_%s_%s_%s_fit_exp_fct_%s_shift_%dperc_%dfact_2.png'
#              % (data_source0, data_source,
#                  percent, time_freq, shift_by_percent,
#                  shift_factor))),
plt.savefig(os.path.join(save_dir,  # filtered_2_neighbors_data
            (r'_test_years.png')),
            papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
