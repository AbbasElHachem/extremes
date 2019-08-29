#!/usr/bin/env python3
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
#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

from scipy.optimize import leastsq

# from _00_additional_functions import (gen_path_df_file)

#==============================================================================
def gen_path_df_file(main_path, start_path_acc, time_freq,
                     data_source, percent, neighbor):
    """ use this funtion to get the path to the dfs for different neighbors"""

    return main_path / (
        r'%sfreq_%s_%s_netatmo_upper_%s_percent_data_considered_neighbor_%d_.csv'
        % (start_path_acc, time_freq, data_source, percent, neighbor))

plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 12})
#==============================================================================
#
#==============================================================================
#%%
#main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

main_dir = Path(r'/run/media/abbas/EL HACHEM/')

main_dir2 = Path(r'/home/abbas/Desktop/data/data/')

data_dir_Netamto_dfs = main_dir2 # / \
#    r'plots_NetAtmo_ppt_DWD_ppt_correlation_'
data_dir_DWD_dfs = main_dir #/ r'plots_DWD_ppt_DWD_ppt_correlation_'

netatmo_path_acc = r'year_allyears_df_comparing_correlations_max_sep_dist_500000_'
dwd_path_Acc = r'year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '95'
time_freq = '60min'

data_source0 = 'DWD'  # reference station 
data_source1 = 'Netatmo'  # target data

data_source = 'dwd'  # compare to station 'netatmo'
#==============================================================================

#==============================================================================
#%%
if data_source1 == 'Netatmo':
    df00 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 0)
    df01 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 1)
    df02 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 2)
    df03 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 3)
    df04 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 4)
    save_dir0 = data_dir_Netamto_dfs
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
#%%
def read_filter_df_corr_return_stns_x_y_vals(df_file, thr_low=0, thr_high=1,
                                             x_col_name='Distance to neighbor',
                                             y_col_name='Bool_Spearman_Correlation'):
    """ read df with correlation values, select all between 0 and 1 and 
        return station_ids, distance_values, correlation_values, and df"""
    in_df = pd.read_csv(df_file, index_col=0, sep=';').dropna(how='any')
    in_df = in_df[(thr_low < in_df.loc[:, y_col_name]) &
                  (in_df.loc[:, y_col_name] < thr_high)]
    stn_ids = in_df.index
    x_vals = in_df.loc[:, x_col_name].values.ravel()
    y_vals = in_df.loc[:, y_col_name].values.ravel()
    return stn_ids, x_vals, y_vals, in_df

#%%

# DWD DATA
stn0, xvals0, yvals0, in_df0 = read_filter_df_corr_return_stns_x_y_vals(df0) 
stn1, xvals1, yvals1, in_df1 = read_filter_df_corr_return_stns_x_y_vals(df1) 
stn2, xvals2, yvals2, in_df2 = read_filter_df_corr_return_stns_x_y_vals(df2) 
stn3, xvals3, yvals3, in_df3 = read_filter_df_corr_return_stns_x_y_vals(df3) 
stn3, xvals4, yvals4, in_df4 = read_filter_df_corr_return_stns_x_y_vals(df4) 

# Netatmo DATA
stn00, xvals00, yvals00, in_df00 = read_filter_df_corr_return_stns_x_y_vals(df00) 
stn01, xvals01, yvals01, in_df01 = read_filter_df_corr_return_stns_x_y_vals(df01) 
stn02, xvals02, yvals02, in_df02 = read_filter_df_corr_return_stns_x_y_vals(df02) 
stn03, xvals03, yvals03, in_df03 = read_filter_df_corr_return_stns_x_y_vals(df03) 
stn03, xvals04, yvals04, in_df04 = read_filter_df_corr_return_stns_x_y_vals(df04) 
# =============================================================================
#%%


#def model(x, a, b, c, d):
#    """ 3degree polynomial function used for fitting and as filter"""
#    return a * x**3 + b * x**2 + c * x + d

# true parameters and a function that takes them
model = lambda x, a, b, c: (a / np.sqrt(x + b) + c)
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
    #print pars, resid.sum() # uncomment to check the iteration parameters
    return np.nan_to_num(resid)



#%%
#def get_x_y_vals_df(df, xcol='Distance to neighbor')
# DWD data
xvals0 = in_df0.loc[:, 'Distance to neighbor'].values.ravel()
yvals0 = in_df0.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

xvals1 = in_df1.loc[:, 'Distance to neighbor'].values.ravel()
yvals1 = in_df1.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

xvals2 = in_df2.loc[:, 'Distance to neighbor'].values.ravel()
yvals2 = in_df2.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

xvals3 = in_df3.loc[:, 'Distance to neighbor'].values.ravel()
yvals3 = in_df3.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

xvals4 = in_df4.loc[:, 'Distance to neighbor'].values.ravel()
yvals4 = in_df4.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

#%%
# Netatmo DATA

#def get_x_y_vals_df(df, xcol='Distance to neighbor')
xvals00 = in_df00.loc[:, 'Distance to neighbor'].values.ravel()
yvals00 = in_df00.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

xvals01 = in_df01.loc[:, 'Distance to neighbor'].values.ravel()
yvals01 = in_df01.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

xvals02 = in_df02.loc[:, 'Distance to neighbor'].values.ravel()
yvals02 = in_df02.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

xvals03 = in_df03.loc[:, 'Distance to neighbor'].values.ravel()
yvals03 = in_df03.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

xvals04 = in_df04.loc[:, 'Distance to neighbor'].values.ravel()
yvals04 = in_df04.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

#%%
# DWD DATA
xall = np.hstack([xvals0, xvals1, xvals2, xvals3, xvals4])
yall = np.hstack([yvals0, yvals1, yvals2, yvals3, yvals4])

# Netatmo DATA
xall0 = np.hstack([xvals00, xvals01, xvals02, xvals03, xvals04])
yall0 = np.hstack([yvals00, yvals01, yvals02, yvals03, yvals04])
#%%

#xall_ = np.insert(arr=xall, obj=0, values=0)
#yall_ = np.insert(arr=yall, obj=0, values=1)
#%%
guesses =[100, 100, 0]

fit_pars, flag = leastsq(func = flipped_resid, x0 = guesses,
                         args = (xall, yall))

y_fit = model(xall, *fit_pars)

# shift curves 5percent downwards
y_fit = y_fit - 0.025*y_fit.min()
#%%
guesses = fit_pars

fit_pars0, flag0 = leastsq(func = flipped_resid, x0 = guesses,
                         args = (xall0, yall0))

y_fit0 = model(xall0, *fit_pars)

# shift curves 5percent downwards
y_fit0 = y_fit0 - 0.025*y_fit0.min()
# %%
plt.ioff()
plt.plot(xall0,
         yall0,
         '.', alpha=0.8, label = 'Test data')
plt.scatter(xall, y_fit, color='r', zorder = 0.9, label = 'Edge')
plt.scatter(xall0, y_fit0, color='g', zorder = 0.9, label = 'Guess')
plt.legend(loc = 'lower left')
plt.show()
#%%

ybelow = yall0[yall0 < y_fit0]