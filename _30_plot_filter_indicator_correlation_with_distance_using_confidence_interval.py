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

from _00_additional_functions import mean_confidence_interval

from _27_plot_indicator_correlation_with_distance_ import gen_path_df_file
from _28_plot_filter_indicator_correlation_with_distance_all_neighbors_combined import (
    func, fit_curve_get_vals_below_curve)


main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

data_dir_Netamto_dfs = main_dir / \
    r'plots_NetAtmo_ppt_DWD_ppt_correlation_'
data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'

netatmo_path_acc = r'year_allyears_df_comparing_correlations_max_sep_dist_500000_'
dwd_path_Acc = r'year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '95'
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


in_df0 = pd.read_csv(df0, index_col=0, sep=';').dropna(how='all')
in_df1 = pd.read_csv(df1, index_col=0, sep=';').dropna(how='all')
in_df2 = pd.read_csv(df2, index_col=0, sep=';').dropna(how='all')
in_df3 = pd.read_csv(df3, index_col=0, sep=';').dropna(how='all')
in_df4 = pd.read_csv(df4, index_col=0, sep=';').dropna(how='all')
# =============================================================================
# TODO: put me in a function
# get all values above 5 percent quantile
in_df0_q = in_df0[in_df0.Bool_Spearman_Correlation.quantile(.01) <
                  in_df0.Bool_Spearman_Correlation]
in_df1_q = in_df1[in_df1.Bool_Spearman_Correlation.quantile(.01) <
                  in_df1.Bool_Spearman_Correlation]
in_df2_q = in_df2[in_df2.Bool_Spearman_Correlation.quantile(.01) <
                  in_df2.Bool_Spearman_Correlation]
in_df3_q = in_df3[in_df3.Bool_Spearman_Correlation.quantile(.01) <
                  in_df3.Bool_Spearman_Correlation]
in_df4_q = in_df4[in_df4.Bool_Spearman_Correlation.quantile(.01) <
                  in_df4.Bool_Spearman_Correlation]
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
# get qunatile values
x0_q = in_df0_q.loc[:, 'Distance to neighbor'].values.ravel()
x1_q = in_df1_q.loc[:, 'Distance to neighbor'].values.ravel()
x2_q = in_df2_q.loc[:, 'Distance to neighbor'].values.ravel()
x3_q = in_df3_q.loc[:, 'Distance to neighbor'].values.ravel()
x4_q = in_df4_q.loc[:, 'Distance to neighbor'].values.ravel()

y0_q = in_df0_q.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y1_q = in_df1_q.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y2_q = in_df2_q.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y3_q = in_df3_q.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y4_q = in_df4_q.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# =============================================================================

# m0, m_h0, mp_h0 = mean_confidence_interval(y0, .95)
#
# m1, m_h1, mp_h1 = mean_confidence_interval(y1, .95)
# m2, m_h2, mp_h2 = mean_confidence_interval(y2, .95)
# m3, m_h3, mp_h3 = mean_confidence_interval(y3, .95)
# m4, m_h4, mp_h4 = mean_confidence_interval(y4, .95)
#
# m_h02 = m_h0 - m_h0 * 0.1
# m_h12 = m_h1 - m_h1 * 0.1
# m_h22 = m_h2 - m_h2 * 0.1
# m_h32 = m_h3 - m_h3 * 0.1
# m_h42 = m_h4 - m_h4 * 0.1
#
# y00 = y0[y0 > m_h0]
# x00 = x0[y0 > m_h0]
#
# y10 = y1[y1 > m_h1]
# x10 = x1[y1 > m_h1]
#
# y20 = y2[y2 > m_h2]
# x20 = x2[y2 > m_h2]
#
# y30 = y3[y3 > m_h3]
# x30 = x3[y3 > m_h3]
#
# y40 = y4[y4 > m_h4]
# x40 = x4[y4 > m_h4]
# # y01 = y0[y0 < mp_h]
# # x01 = x0[y0 < mp_h]
#
#
# y002 = y0[y0 > m_h02]
# x002 = x0[y0 > m_h02]
#
# y102 = y1[y1 > m_h12]
# x102 = x1[y1 > m_h12]
#
# y202 = y2[y2 > m_h22]
# x202 = x2[y2 > m_h22]
#
# y302 = y3[y3 > m_h32]
# x302 = x3[y3 > m_h32]
#
# y402 = y4[y4 > m_h42]
# x402 = x4[y4 > m_h42]
#
#
# plt.hist(y0, bins=15, color='blue', edgecolor='black')
# plt.hist(y00, bins=15, color='red', edgecolor='black')
# plt.hist(y01, bins=15, color='g', edgecolor='black')

# plt.scatter(x0, y0, c='r', alpha=0.5)

# plt.show()
# %% apply filter for every neighbor seperatly

(y_fitted_shifted0, xvals_below_curve0,
 yvals_below_curve0, xvals_above_curve0,
 yvals_above_curve0, stnbelow0, stnabove0) = fit_curve_get_vals_below_curve(
    x0, y0, func, s0)

(y_fitted_shifted1, xvals_below_curve1,
 yvals_below_curve1, xvals_above_curve1,
 yvals_above_curve1, stnbelow1, stnabove1) = fit_curve_get_vals_below_curve(
    x1, y1, func, s1)
(y_fitted_shifted2, xvals_below_curve2,
 yvals_below_curve2, xvals_above_curve2,
 yvals_above_curve2, stnbelow2, stnabove2) = fit_curve_get_vals_below_curve(
    x2, y2, func, s2)

(y_fitted_shifted3, xvals_below_curve3,
 yvals_below_curve3, xvals_above_curve3,
 yvals_above_curve3, stnbelow3, stnabove3) = fit_curve_get_vals_below_curve(
    x3, y3, func, s3)

(y_fitted_shifted4, xvals_below_curve4,
 yvals_below_curve4, xvals_above_curve4,
 yvals_above_curve4, stnbelow4, stnabove4) = fit_curve_get_vals_below_curve(
    x4, y4, func, s4)

# =============================================================================

# TODO: CHECK AGAIN intersect all stations

stns_keep_all = np.intersect1d(np.intersect1d(np.intersect1d(np.intersect1d(
    stnabove0, stnabove1),
    stnabove2), stnabove3), stnabove4)

stns_keep_al_sr = pd.Series(stns_keep_all)
stns_keep_al_sr.to_csv(
    (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
        r'\filter_Netamo_data_basedon_indicator_correlation'
        r'\keep_stns_all_neighbors_combined_%s_per_.csv' % percent),
    sep=';')
# %%
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
#
#==============================================================================
plt.ioff()
plt.figure(figsize=(16, 12), dpi=300)

plt.scatter(x0, y0, c='r', alpha=0.5)
plt.scatter(x1, y1, c='r', alpha=0.5)
plt.scatter(x2, y2, c='r', alpha=0.5)
plt.scatter(x3, y3, c='r', alpha=0.5)
plt.scatter(x4, y4, c='r', alpha=0.5)

# plt.scatter(x00, y00, c='c', alpha=0.5)
# plt.show()
# plt.scatter(x10, y10, c='c', alpha=0.5)
# plt.scatter(x20, y20, c='c', alpha=0.5)
# plt.scatter(x30, y30, c='c', alpha=0.5)
# plt.scatter(x40, y40, c='c', alpha=0.5)

plt.scatter(x0_q, y0_q, c='k', alpha=0.75)
plt.scatter(x1_q, y1_q, c='k', alpha=0.75)
plt.scatter(x2_q, y2_q, c='k', alpha=0.75)
plt.scatter(x3_q, y3_q, c='k', alpha=0.75)
plt.scatter(x4_q, y4_q, c='k', alpha=0.75)

plt.scatter(x0_abv, y0_abv, c='r', alpha=0.5,
            marker='x', label='First Neighbor', s=28)
plt.scatter(x1_abv, y1_abv, c='b', alpha=0.5,
            marker='.', label='Second Neighbor', s=28)
plt.scatter(x2_abv, y2_abv, c='g', alpha=0.5,
            marker='d', label='Third Neighbor', s=28)
plt.scatter(x3_abv, y3_abv, c='orange', alpha=0.5,
            marker='1', label='Fourth Neighbor', s=28)
plt.scatter(x4_abv, y4_abv, c='m', alpha=0.5,
            marker='1', label='Fifth Neighbor', s=28)


# plt.scatter(x002, y002, c='k', alpha=0.5)
# plt.scatter(x102, y102, c='k', alpha=0.5)
# plt.scatter(x202, y202, c='k', alpha=0.5)
# plt.scatter(x302, y302, c='k', alpha=0.5)
# plt.scatter(x402, y402, c='k', alpha=0.5)
# plt.scatter(xvals_below_curve0, yvals_below_curve0, c='k', alpha=0.95,
#             marker='x', label='First Neighbor', s=28)
# plt.scatter(xvals_below_curve1, yvals_below_curve1, c='k', alpha=0.95,
#             marker='.', label='Second Neighbor', s=28)
# plt.scatter(xvals_below_curve2, yvals_below_curve2, c='k', alpha=0.95,
#             marker='d', label='Third Neighbor', s=28)
# plt.scatter(xvals_below_curve3, yvals_below_curve3, c='k', alpha=0.95,
#             marker='1', label='Fourth Neighbor', s=28)
# plt.scatter(xvals_below_curve4, yvals_below_curve4, c='k', alpha=0.95,
#             marker='1', label='Fifth Neighbor', s=28)

# plt.scatter(x0, y_fitted_shifted0, c='darkred', alpha=0.5,
#             marker='x', label='Fitted curve 1', s=20)
#
# plt.scatter(x1, y_fitted_shifted1, c='darkblue', alpha=0.5,
#             marker='x', label='Fitted curve 2', s=20)
#
# plt.scatter(x2, y_fitted_shifted2, c='darkgreen', alpha=0.5,
#             marker='x', label='Fitted curve 3', s=20)
#
# plt.scatter(x3, y_fitted_shifted3, c='darkorange', alpha=0.5,
#             marker='x', label='Fitted curve 4', s=20)
# plt.scatter(x4, y_fitted_shifted4, c='purple', alpha=0.5,
#             marker='x', label='Fitted curve 5', s=20)


plt.xlim([0, max([x3_q.max(), x4_q.max()]) + 1000])
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
# plt.tight_layout()
plt.title('Keeping %s %d stations: Indicator correlation with distance'
          ' for upper %s percent of data values'
          % (data_source0, stns_keep_all.shape[0], percent))
plt.savefig(save_dir /
            (r'_%s_%s_%s_percent_indic_corr_freq_%s_each_neighbors7.png'
             % (data_source0, data_source, percent, time_freq)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
