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

This is used as a filter to remove the 'bad' Netamo staions 


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

from pathlib import Path

from _00_additional_functions import (
    gen_path_df_file, read_filter_df_corr_return_stns_x_y_vals)

plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 14})
#==============================================================================
#
#==============================================================================
main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_DWD_ppt_correlation_'

assert data_dir_Netamto_dfs.exists(), 'Wrong Netatmo path'


data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'
assert data_dir_DWD_dfs.exists(), 'Wrong dwd path'

data_dir_Netamto_netatmo_dfs = main_dir / \
    r'plots_NetAtmo_ppt_Netatmo_ppt_correlation_'

assert data_dir_Netamto_netatmo_dfs.exists(), 'Wrong Netatmo Netatmo path'


# allyears pearson_  new_method new_method_pearson_
netatmo_path_acc = r'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_30000_'
dwd_path_Acc = r'pearson_year_allyears_df_dwd_correlations'
# freq_60min_dwd_netatmo_upper_99_0_percent_data_considered_neighbor_0_

path_to_netatmo_gd_stns_file = data_dir_Netamto_dfs / \
    r'keep_stns_all_neighbor_99_per_60min_s0_comb.csv'

#assert path_to_netatmo_gd_stns_file.exists(), 'wrong netatmo good stns file'

# def percentage threshold, time frequency and data source
percent = '99_0'
time_freq = '60min'  # '720min', 1440min, '480min', '360min', '180min', '120min'
# '60min'
data_source0 = 'Netatmo'  # 'DWD'  # 'Netatmo'  #   # reference station 'Netatmo'
data_source = 'dwd'  # 'dwd'  # 'netatmo'  #   # compare to station 'netatmo'

use_good_netatmo_stns = False
use_filtered_data = False
filtered_percent = '99_0'

save_acc = ''
# =============================================================================


if data_source0 == 'Netatmo' and data_source == 'dwd':
    # for Netatmo stations neighbors start from 0 !
    df0 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 0,
                           use_filtered_data=use_filtered_data,
                           filter_percent=filtered_percent)
    df1 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 1,
                           use_filtered_data=use_filtered_data,
                           filter_percent=filtered_percent)
    df2 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 2,
                           use_filtered_data=use_filtered_data,
                           filter_percent=filtered_percent)
#     df3 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 3, use_filtered_data)
#     df4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 4, use_filtered_data)
#     df5 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 5, use_filtered_data)
#     df6 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 6, use_filtered_data)
#     df7 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
#                            data_source, percent, 7, use_filtered_data)
    save_dir = data_dir_Netamto_dfs

data_source0 = 'DWD'
percent = '99'
if data_source0 == 'DWD' and data_source == 'dwd':
    # for DWD stations neighbors start from 1 not 0  (1 is first)!
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

if data_source0 == 'Netatmo' and data_source == 'netatmo':
    # for Netatmo stations neighbors start from 0 !
    df0 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
                           netatmo_path_acc, time_freq,
                           data_source, percent, 0,
                           use_filtered_data=use_filtered_data,
                           filter_percent=filtered_percent)
    df1 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
                           netatmo_path_acc, time_freq,
                           data_source, percent, 1,
                           use_filtered_data=use_filtered_data,
                           filter_percent=filtered_percent)
#     df2 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
#                            netatmo_path_acc, time_freq,
#                            data_source, percent, 2)
#     df3 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
#                            netatmo_path_acc, time_freq,
#                            data_source, percent, 3)
#     df4 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
#                            netatmo_path_acc, time_freq,
#                            data_source, percent, 4)
    save_dir = data_dir_Netamto_netatmo_dfs
# =============================================================================

# for DWD stations neighbors start from 1 not 0  (1 is first)!


in_df0_dwd = pd.read_csv(df0_dwd, index_col=0, sep=';').dropna(how='all')
in_df1_dwd = pd.read_csv(df1_dwd, index_col=0, sep=';').dropna(how='all')
# in_df2_dwd = pd.read_csv(df2_dwd, index_col=0, sep=';').dropna(how='all')
# Netatmo

#'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_30000_freq_60min_dwd_netatmo_upper_99_0_percent_data_considered_neighbor_0_'
in_df0 = pd.read_csv(df0, index_col=0, sep=';').dropna(how='all')
in_df1 = pd.read_csv(df1, index_col=0, sep=';').dropna(how='all')
# in_df2 = pd.read_csv(df2, index_col=0, sep=';').dropna(how='all')
# in_df3 = pd.read_csv(df3, index_col=0, sep=';').dropna(how='all')
# in_df4 = pd.read_csv(df4, index_col=0, sep=';').dropna(how='all')
# in_df5 = pd.read_csv(df5, index_col=0, sep=';').dropna(how='all')
# in_df6 = pd.read_csv(df6, index_col=0, sep=';').dropna(how='all')
# in_df7 = pd.read_csv(df7, index_col=0, sep=';').dropna(how='all')

if use_good_netatmo_stns:
    df_good_stns = pd.read_csv(path_to_netatmo_gd_stns_file, sep=';',
                               index_col=0).dropna(how='any')
    # assert np.all(df_good_stns.values.ravel()) in in_df0.index
    in_df0 = in_df0.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df1 = in_df1.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df2 = in_df2.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df3 = in_df3.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df4 = in_df4.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df5 = in_df2.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df6 = in_df3.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df7 = in_df4.loc[df_good_stns.values.ravel(), :].dropna(how='all')
    save_acc = 'filtered_using_stns_abv_curve_all_neighbor'
# =============================================================================


x0_dwd0 = in_df0_dwd.loc[:, 'Distance to neighbor'].values.ravel()
x1_dwd0 = in_df1_dwd.loc[:, 'Distance to neighbor'].values.ravel()
# x2_dwd0 = in_df2_dwd.loc[:, 'Distance to neighbor'].values.ravel()

y0_dwd0 = in_df0_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y1_dwd0 = in_df1_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y2_dwd0 = in_df2_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()


in_df0 = in_df0[in_df0['Bool_Pearson_Correlation_Netatmo_DWD'] > 0.2]
in_df1 = in_df1[in_df1['Bool_Pearson_Correlation_Netatmo_DWD'] > y1_dwd0.min()]
in_df0 = in_df0[in_df0['Bool_Pearson_Correlation_DWD_DWD'] > 0.2]
in_df1 = in_df1[in_df1['Bool_Pearson_Correlation_DWD_DWD'] > y1_dwd0.min()]

# in_df0_dwd = in_df0_dwd[in_df0_dwd['Bool_Spearman_Correlation'] > 0.2]
# in_df1_dwd = in_df1_dwd[in_df1_dwd['Bool_Spearman_Correlation'] > 0.2]

stns_keep_all_final_new = in_df0.index.intersection(in_df1.index)

# stns_keep_all_final_new = in_df0.index

#cmn_stns = in_df0.index.intersection(in_df1.index)
x0 = in_df0.loc[stns_keep_all_final_new, 'Distance to neighbor'].values.ravel()
x1 = in_df1.loc[stns_keep_all_final_new, 'Distance to neighbor'].values.ravel()
# x2 = in_df2.loc[:, 'Distance to neighbor'].values.ravel()
# x3 = in_df3.loc[:, 'Distance to neighbor'].values.ravel()
# x4 = in_df4.loc[:, 'Distance to neighbor'].values.ravel()
# x5 = in_df5.loc[:, 'Distance to neighbor'].values.ravel()
# x6 = in_df6.loc[:, 'Distance to neighbor'].values.ravel()
# x7 = in_df7.loc[:, 'Distance to neighbor'].values.ravel()

y0 = in_df0.loc[stns_keep_all_final_new,
                'Bool_Pearson_Correlation_Netatmo_DWD'].values.ravel()
y1 = in_df1.loc[stns_keep_all_final_new,
                'Bool_Pearson_Correlation_Netatmo_DWD'].values.ravel()

# y0_dwd = in_df0.loc[:,
#                     'Bool_Pearson_Correlation_DWD_DWD'].values.ravel()
# y1_dwd = in_df1.loc[:,
#                     'Bool_Pearson_Correlation_DWD_DWD'].values.ravel()

# y0 = in_df0.loc[:,
#                 'Bool_Spearman_Correlation'].values.ravel()
# y1 = in_df1.loc[:,
#                 'Bool_Spearman_Correlation'].values.ravel()

# y2 = in_df2.loc[:,
#                 'Bool_Spearman_Correlation'].values.ravel()
# y0_sens = in_df0.loc[stns_keep_all_final_new,
#                      'Bool_Sensored_Pearson_Correlation'].values.ravel()
# y1_sens = in_df1.loc[stns_keep_all_final_new,
#                      'Bool_Sensored_Pearson_Correlation'].values.ravel()

# y2 = in_df2.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y3 = in_df3.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y4 = in_df4.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y5 = in_df5.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y6 = in_df6.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y7 = in_df7.loc[:, 'Bool_Spearman_Correlation'].values.ravel()


# stns_keep_al_sr = pd.DataFrame(data=stns_keep_all_final_new,
#                                columns=['Stations'])
#
# stns_keep_al_sr.to_csv(
#     (save_dir /
#         (r'keep_stns_all_neighbor_%s_per_%s_s0_1st.csv'
#          % (percent, time_freq))),
#     sep=';')


# s0, x0, y0, in_df0 = read_filter_df_corr_return_stns_x_y_vals(df0)
# =============================================================================

plt.ioff()
plt.figure(figsize=(12, 8), dpi=300)

plt.scatter(x0, y0, c='r', alpha=0.5, marker='x',
            label='First Neighbor %d Pairs' % y0.shape[0], s=34)

plt.scatter(x0_dwd0, y0_dwd0, c='orange', alpha=0.5, marker='o',
            label='DWD First Neighbor', s=34)

plt.scatter(x1, y1, c='blue', alpha=0.5, marker='.',
            label='Second Neighbor %d Pairs' % y1.shape[0], s=34)


# plt.scatter(x0, y0_dwd, c='orange', alpha=0.35, marker='.',
#             label='First Neighbor %d Pairs' % y0.shape[0], s=24)
#
# # plt.scatter(x0, y0_dwd, c='orange', alpha=0.5, marker='o',
# #             label='DWD First Neighbor', s=34)
#
plt.scatter(x1_dwd0, y1_dwd0, c='g', alpha=0.35, marker='1',
            label='DWD Second Neighbor', s=34)

# plt.scatter(x2, y2, c='green', alpha=0.35, marker='d',
#             label='Third Neighbor Stn nbr %d' % y2.shape[0], s=34)

# plt.scatter(x1, y1_dwd, c='g', alpha=0.5, marker='*',
#             label=' DWD Second Neighbor', s=34)

# plt.scatter(x0_dwd0, y0_dwd0, c='k', alpha=0.75, marker='.',
#             label='DWD First Neighbor', s=34)
# #
# plt.scatter(x1_dwd0, y1_dwd0, c='k', alpha=0.75, marker='.',
#             label='DWD Second Neighbor', s=34)

# plt.scatter(x2_dwd0, y2_dwd0, c='k', alpha=0.95, marker='d',
#             label=' DWD Third Neighbor', s=34)


# plt.scatter(x3, y3, c='darkorange', alpha=0.5, marker='*',
#             label='Fourth Neighbor Stn nbr %d' % y3.shape[0], s=34)
# plt.scatter(x4, y4, c='m', alpha=0.5, marker='+',
#             label='Fifth Neighbor Stn nbr %d' % y4.shape[0], s=34)
# plt.scatter(x5, y5, c='darkgreen', alpha=0.5, marker='1',
#             label='Sixth Neighbor Stn nbr %d' % y5.shape[0], s=34)
# plt.scatter(x6, y6, c='darkblue', alpha=0.5, marker='X',
#             label='Seventh Neighbor Stn nbr %d' % y6.shape[0], s=34)
# plt.scatter(x7, y7, c='c', alpha=0.5, marker='8',
#             label='Eighth Neighbor Stn nbr %d' % y7.shape[0], s=34)
plt.xlim([0, max(x0.max(), 25000) + 100])
plt.xticks(np.arange(0, 25000 + 100, 5000))  # x1.max()
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
plt.tight_layout()

if data_source == 'dwd':
    data_source = 'DWD'
if use_filtered_data:
    save_acc = save_acc + 'filtered_' + str(filtered_percent)

# if temp_freq == '60min':
#     temp_freq = 'hourly'

# plt.title('Indicator correlation for %s %s stations for %s precipitation, \n '
#          ' for upper %s percent of data values '
#          % (data_source0, data_source, time_freq, percent))
plt.savefig(save_dir /
            (r'new_pearson__%s_%s_%s_percent_indic_corr_freq_%s_%s_10_2.png'
             % (data_source0, data_source, percent, time_freq, save_acc)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
plt.show()
