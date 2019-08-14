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

from scipy.optimize import curve_fit
from pathlib import Path

from _24_plot_indicator_correlation_with_distance_ import gen_path_df_file

main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_NetAtmo_temperature'
data_dir_DWD_dfs = main_dir / r'plots_indicator_correlations_DWD_DWD_ppt_'

netatmo_path_acc = r'year_allyears_df_comparing_correlations_max_sep_dist_500000_'
dwd_path_Acc = r'year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '95'
time_freq = '60min'

data_source0 = 'DWD'  # reference station 'Netatmo'
data_source = 'dwd'  # compare to station 'netatmo'
#==============================================================================
#
#==============================================================================

# %%
# def func(x, a, b, c):
#    return a * np.exp(-b * x) + c


def func(x, a, b, c, d):
    ''' 3degree polynomial function used for fitting and as filter'''
    return a * x**3 + b * x**2 + c * x + d


def fit_curve_get_vals_below_curve(x, y, func, stns):
    ''' fit function to data and shifted 10% downwards'''

    popt, _ = curve_fit(func, x, y)
    y_fitted = func(x, *popt)

    lower_bound = y_fitted.max() * 0.10
    y_fitted_shifted = y_fitted - lower_bound

    xvals_below_curve = x[np.where(y <= y_fitted_shifted)]
    yvals_below_curve = y[np.where(y <= y_fitted_shifted)]

    stns_below_curve = stns[np.where(y <= y_fitted_shifted)]

    xvals_above_curve = x[np.where(y > y_fitted_shifted)]
    yvals_above_curve = y[np.where(y > y_fitted_shifted)]

    stns_above_curve = stns[np.where(y > y_fitted_shifted)]

    return (y_fitted_shifted, xvals_below_curve, yvals_below_curve,
            xvals_above_curve, yvals_above_curve, stns_below_curve,
            stns_above_curve)


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


in_df0 = pd.read_csv(df0, index_col=0, sep=';').dropna(how='all')
in_df1 = pd.read_csv(df1, index_col=0, sep=';').dropna(how='all')
in_df2 = pd.read_csv(df2, index_col=0, sep=';').dropna(how='all')
in_df3 = pd.read_csv(df3, index_col=0, sep=';').dropna(how='all')
in_df4 = pd.read_csv(df4, index_col=0, sep=';').dropna(how='all')

# =============================================================================


s0 = in_df0.index
s1 = in_df1.index
s2 = in_df2.index
s3 = in_df3.index
s4 = in_df4.index

x0 = in_df0.loc[:, 'Distance to neighbor'].values.ravel()
x1 = in_df1.loc[:, 'Distance to neighbor'].values.ravel()
x2 = in_df2.loc[:, 'Distance to neighbor'].values.ravel()
x3 = in_df3.loc[:, 'Distance to neighbor'].values.ravel()
x4 = in_df4.loc[:, 'Distance to neighbor'].values.ravel()

y0 = in_df0.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y1 = in_df1.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y2 = in_df2.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y3 = in_df3.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y4 = in_df4.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

# =============================================================================

# combined all stns, distances and correlation data into arrays
stns = np.concatenate((s0, s1, s2, s3, s4))
xs = np.concatenate((x0, x1, x2, x3, x4))
ys = np.concatenate((y0, y1, y2, y3, y4))

# remove nans
stns = stns[np.where(ys >= 0)]
xs = xs[np.where(ys >= 0)]
ys = ys[ys >= 0]

# fit function to all data combined
(yfs, xvals_below_curves, yvals_below_curves,
 xvals_above_curves, yvals_above_curves,
 stns_below_curves, stns_above_curves) = fit_curve_get_vals_below_curve(
    xs, ys, func, stns)

# =============================================================================

stns_to_keep = np.unique(stns_above_curves)
stns_to_remove = np.unique(stns_below_curves)

# TODO: NEED FIXING ASK PROF
stns_keep_all = np.intersect1d(
    np.intersect1d(np.intersect1d(np.intersect1d(np.intersect1d(
        stns_to_keep, s0),
        s1), s2), s3), s4)

# stns_keep_al_sr = pd.Series(stns_keep_all)
# stns_keep_al_sr.to_csv(
#     (main_dir /
#         r'filter_Netamo_data_basedon_indicator_correlation' /
#         r'keep_stns_all_neighbors_combined_%s_per_.csv' % percent),
#     sep=';')
# =============================================================================

x0_abv = in_df0.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y0_abv = in_df0.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

x1_abv = in_df1.loc[stns_keep_all, 'Distance to neighbor']  # .dropna().values
y1_abv = in_df1.loc[stns_keep_all,
                    'Bool_Spearman_Correlation']  # .dropna().values

x2_abv = in_df2.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y2_abv = in_df2.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

x3_abv = in_df3.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y3_abv = in_df3.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

x4_abv = in_df4.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y4_abv = in_df4.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

# TODO: Add Additional filter on stations
# =============================================================================

plt.ioff()
plt.figure(figsize=(12, 8), dpi=100)


plt.scatter(x0, y0, c='r', alpha=0.5,
            marker='x', label='First Neighbor', s=28)
plt.scatter(x1, y1, c='orange', alpha=0.5,
            marker='.', label='Second Neighbor', s=28)
plt.scatter(x2, y2, c='g', alpha=0.5,
            marker='d', label='Third Neighbor', s=28)
plt.scatter(x3, y3, c='b', alpha=0.5,
            marker='1', label='Fourth Neighbor', s=28)
plt.scatter(x4, y4, c='b', alpha=0.5,
            marker='1', label='Fourth Neighbor', s=28)

plt.scatter(xs, yfs, label="Polynomial function",
            alpha=.25, c='m', marker='*', s=20)
plt.scatter(xvals_above_curves, yvals_above_curves, s=20,
            label="Stns above Curve%d" % stns_to_keep.shape[0],
            alpha=.25, c='grey', marker='d')

plt.scatter(xvals_below_curves, yvals_below_curves, s=20,
            label="Below Curve %d" % stns_to_remove.shape[0],
            alpha=.25, c='k', marker='X')

#plt.scatter([],[], label='Number of stations %d' % dfs_all.shape[0], marker='.')
plt.xlim([0, x4.max() + 1000])
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
plt.tight_layout()
# plt.axis('equal')
plt.title('Keeping %s %d stations: Indicator correlation with distance'
          ' for upper %s percent of data values'
          % (data_source0, stns_to_keep.shape[0], percent))
plt.savefig(save_dir /
            (r'_%s_%s_%s_percent_indic_corr_freq_%s_all_neighbors.png'
             % (data_source0, data_source, percent, time_freq)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
# plt.show()
plt.close()
