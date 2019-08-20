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

from _00_additional_functions import (gen_path_df_file)


plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 12})
#==============================================================================
#
#==============================================================================
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


in_df0 = pd.read_csv(df0, index_col=0, sep=';').dropna(how='any')
in_df1 = pd.read_csv(df1, index_col=0, sep=';').dropna(how='any')
in_df2 = pd.read_csv(df2, index_col=0, sep=';').dropna(how='any')
in_df3 = pd.read_csv(df3, index_col=0, sep=';').dropna(how='any')
in_df4 = pd.read_csv(df4, index_col=0, sep=';').dropna(how='any')

in_df0 = in_df0[in_df0.Bool_Spearman_Correlation < 1]
in_df1 = in_df1[in_df1.Bool_Spearman_Correlation < 1]
in_df2 = in_df2[in_df2.Bool_Spearman_Correlation < 1]
in_df3 = in_df3[in_df3.Bool_Spearman_Correlation < 1]
in_df4 = in_df4[in_df4.Bool_Spearman_Correlation < 1]

# =============================================================================
# TODO: put me in a function
# get all values above 5 percent quantile
in_df0_q = in_df0[in_df0.Bool_Spearman_Correlation.quantile(.2) <
                  in_df0.Bool_Spearman_Correlation]
in_df1_q = in_df1[in_df1.Bool_Spearman_Correlation.quantile(.25) <
                  in_df1.Bool_Spearman_Correlation]
in_df2_q = in_df2[in_df2.Bool_Spearman_Correlation.quantile(.25) <
                  in_df2.Bool_Spearman_Correlation]
in_df3_q = in_df3[in_df3.Bool_Spearman_Correlation.quantile(.25) <
                  in_df3.Bool_Spearman_Correlation]
in_df4_q = in_df4[in_df4.Bool_Spearman_Correlation.quantile(.25) <
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
# #
# # m_h02 = m_h0 - m_h0 * 0.1
# # m_h12 = m_h1 - m_h1 * 0.1
# # m_h22 = m_h2 - m_h2 * 0.1
# # m_h32 = m_h3 - m_h3 * 0.1
# # m_h42 = m_h4 - m_h4 * 0.1
# #
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
# y01 = y0[y0 < mp_h]
# x01 = x0[y0 < mp_h]
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


#==============================================================================
#
#==============================================================================
plt.ioff()
plt.figure(figsize=(16, 12), dpi=300)

# plt.scatter(x00, y00, c='r', alpha=0.5,
#             marker='x', label='First Neighbor', s=28)
# plt.scatter(x10, y10, c='b', alpha=0.5,
#             marker='.', label='Second Neighbor', s=28)
# plt.scatter(x20, y20, c='g', alpha=0.5,
#             marker='d', label='Third Neighbor', s=28)
# plt.scatter(x30, y30, c='darkorange', alpha=0.5,
#             marker='1', label='Fourth Neighbor', s=28)
# plt.scatter(x40, y40, c='m', alpha=0.5,
#             marker='1', label='Fifth Neighbor', s=28)

plt.scatter(x0_q, y0_q, c='r', alpha=0.5,
            marker='x', label='First Neighbor', s=28)
plt.scatter(x1_q, y1_q, c='b', alpha=0.5,
            marker='.', label='Second Neighbor', s=28)
plt.scatter(x2_q, y2_q, c='g', alpha=0.5,
            marker='d', label='Third Neighbor', s=28)
plt.scatter(x3_q, y3_q, c='darkorange', alpha=0.5,
            marker='1', label='Fourth Neighbor', s=28)
plt.scatter(x4_q, y4_q, c='m', alpha=0.5,
            marker='1', label='Fifth Neighbor', s=28)


plt.xlim([0, max([x3.max(), x4.max()]) + 1000])
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
# plt.tight_layout()
plt.title('Keeping %s %d stations: Indicator correlation with distance'
          ' for upper %s percent of data values'
          % (data_source0, x2_q.shape[0], percent))

plt.savefig(save_dir /
            (r'_%s_%s_%s_percent_indic_corr_freq_%s_filtered_95per_conf_int_.png'
             % (data_source0, data_source, percent, time_freq)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
