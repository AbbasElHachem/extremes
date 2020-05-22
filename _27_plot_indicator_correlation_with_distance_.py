# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Plot Pearson Indicator Correlation with Distance
Purpose: Filter PWS data based on Pearson indicator correlation

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

plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'axes.labelsize': 26})
#==============================================================================
#
#==============================================================================
# for BW
#main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

#data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_DWD_ppt_correlation_'
# r'plots_NetAtmo_1Deg_ppt_DWD_ppt_correlation_'

# for RLP
data_dir_Netamto_dfs = Path(
    r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\indicator_correlation')

# data_dir_Netamto_dfs = Path(
# r'/run/media/abbas/EL Hachem
# 2019/home_office/2020_10_03_Rheinland_Pfalz/indicator_correlation')

data_dir_DWD_dfs = data_dir_Netamto_dfs

assert data_dir_Netamto_dfs.exists(), 'Wrong Netatmo path'


# data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'
# assert data_dir_DWD_dfs.exists(), 'Wrong dwd path'


# data_dir_Netamto_netatmo_dfs = main_dir / \
#     r'plots_NetAtmo_ppt_Netatmo_ppt_correlation_'

# assert data_dir_Netamto_netatmo_dfs.exists(), 'Wrong Netatmo Netatmo path'


# allyears pearson_  new_method new_method_pearson_
netatmo_path_acc_b4 = r'prs_allyears_df_corr_sep_dist_1000000_'

# netatmo_path_acc_after = r'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_1000000_'

netatmo_path_acc_after = r'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_1000000_'

dwd_path_Acc = r'pearson_year_allyears_df_dwd_correlations'
# freq_60min_dwd_netatmo_upper_99_0_percent_data_considered_neighbor_0_
#
# path_to_netatmo_gd_stns_file = data_dir_Netamto_dfs / \
#     r'keep_stns_all_neighbor_99_per_60min_s0_1st.csv'

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


# if data_source0 == 'Netatmo' and data_source == 'dwd':
# for Netatmo stations neighbors start from 0 !
df0_b4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_b4, time_freq,
                          data_source, percent, 0,
                          use_filtered_data=use_filtered_data,
                          filter_percent=filtered_percent)
# df1_b4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_b4, time_freq,
#                           data_source, percent, 1,
#                           use_filtered_data=use_filtered_data,
#                           filter_percent=filtered_percent)
# df2_b4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_b4, time_freq,
#                           data_source, percent, 2,
#                           use_filtered_data=use_filtered_data,
#                           filter_percent=filtered_percent)

# after filtering
df0_after = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_after, time_freq,
                             data_source, percent, 0,
                             use_filtered_data=use_filtered_data,
                             filter_percent=filtered_percent)
# df1_after = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_after, time_freq,
#                              data_source, percent, 1,
#                              use_filtered_data=use_filtered_data,
#                              filter_percent=filtered_percent)
# df2_after = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc_after, time_freq,
#                              data_source, percent, 2,
#                              use_filtered_data=use_filtered_data,
#                              filter_percent=filtered_percent)
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


# percent = '95_0'

# if data_source0 == 'DWD' and data_source == 'dwd':
# for DWD stations neighbors start from 1 not 0  (1 is first)!
df0_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 1)
df1_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 2)
# df2_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 3)


# df3_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 4)
# df4_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 5)
# df5_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 6)
# df6_dwd = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
#                            data_source, percent, 7)
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


# if data_source0 == 'Netatmo' and data_source == 'netatmo':
#     # for Netatmo stations neighbors start from 0 !
#     df0 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
#                            netatmo_path_acc, time_freq,
#                            data_source, percent, 0,
#                            use_filtered_data=use_filtered_data,
#                            filter_percent=filtered_percent)
#     df1 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
#                            netatmo_path_acc, time_freq,
#                            data_source, percent, 1,
#                            use_filtered_data=use_filtered_data,
#                            filter_percent=filtered_percent)
#     df2 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
#                            netatmo_path_acc, time_freq,
#                            data_source, percent, 2)
#     df3 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
#                            netatmo_path_acc, time_freq,
#                            data_source, percent, 3)
#     df4 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
#                            netatmo_path_acc, time_freq,
#                            data_source, percent, 4)
#     save_dir = data_dir_Netamto_netatmo_dfs

save_dir = data_dir_Netamto_dfs
# =============================================================================

# for DWD stations neighbors start from 1 not 0  (1 is first)!


in_df0_dwd = pd.read_csv(df0_dwd, index_col=0, sep=';').dropna(how='all')
in_df1_dwd = pd.read_csv(df1_dwd, index_col=0, sep=';').dropna(how='all')
# in_df2_dwd = pd.read_csv(df2_dwd, index_col=0, sep=';').dropna(how='all')
# in_df3_dwd = pd.read_csv(df3_dwd, index_col=0, sep=';').dropna(how='all')
# in_df4_dwd = pd.read_csv(df4_dwd, index_col=0, sep=';').dropna(how='all')
# in_df5_dwd = pd.read_csv(df5_dwd, index_col=0, sep=';').dropna(how='all')
# in_df6_dwd = pd.read_csv(df6_dwd, index_col=0, sep=';').dropna(how='all')
# Netatmo

#'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_30000_freq_60min_dwd_netatmo_upper_99_0_percent_data_considered_neighbor_0_'
in_df0_b4 = pd.read_csv(df0_b4, index_col=0, sep=';').dropna(how='all')
# in_df1_b4 = pd.read_csv(df1_b4, index_col=0, sep=';').dropna(how='all')
# in_df2_b4 = pd.read_csv(df2_b4, index_col=0, sep=';').dropna(how='all')

in_df0_after = pd.read_csv(df0_after, index_col=0, sep=';').dropna(how='all')
# in_df1_after = pd.read_csv(df1_after, index_col=0, sep=';').dropna(how='all')
# in_df2_after = pd.read_csv(df2_after, index_col=0, sep=';').dropna(how='all')

# if use_good_netatmo_stns:
#     df_good_stns = pd.read_csv(path_to_netatmo_gd_stns_file, sep=';',
#                                index_col=0).dropna(how='any')
#     # assert np.all(df_good_stns.values.ravel()) in in_df0.index
#     in_df0 = in_df0.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df1 = in_df1.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df2 = in_df2.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df3 = in_df3.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df4 = in_df4.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df5 = in_df2.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df6 = in_df3.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     in_df7 = in_df4.loc[df_good_stns.values.ravel(), :].dropna(how='all')
#     save_acc = 'filtered_using_stns_abv_curve_all_neighbor'
# =============================================================================


x0_dwd0 = in_df0_dwd.loc[:, 'Distance to neighbor'].values.ravel()
x1_dwd0 = in_df1_dwd.loc[:, 'Distance to neighbor'].values.ravel()
# x2_dwd0 = in_df2_dwd.loc[:, 'Distance to neighbor'].values.ravel()

# x3_dwd0 = in_df3_dwd.loc[:, 'Distance to neighbor'].values.ravel()
# x4_dwd0 = in_df4_dwd.loc[:, 'Distance to neighbor'].values.ravel()
# x5_dwd0 = in_df5_dwd.loc[:, 'Distance to neighbor'].values.ravel()


y0_dwd0 = in_df0_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y1_dwd0 = in_df1_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y2_dwd0 = in_df2_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

# y3_dwd0 = in_df3_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y4_dwd0 = in_df4_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
# y5_dwd0 = in_df5_dwd.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

x0_b4 = in_df0_b4.loc[
    :, 'Distance to neighbor'].values.ravel()
y0_b4 = in_df0_b4.loc[:,
                      'Bool_Spearman_Correlation'].values.ravel()
#
# x1_b4 = in_df1_b4.loc[
#     :, 'Distance to neighbor'].values.ravel()
# y1_b4 = in_df1_b4.loc[:,
#                       'Bool_Spearman_Correlation'].values.ravel()
#
# x2_b4 = in_df2_b4.loc[
#     :, 'Distance to neighbor'].values.ravel()
# y2_b4 = in_df2_b4.loc[:,
#                       'Bool_Spearman_Correlation'].values.ravel()

in_df0_after = in_df0_after[
    in_df0_after['Bool_Pearson_Correlation_Netatmo_DWD'] > 0.0]
# in_df1_after = in_df1_after[
#     in_df1_after['Bool_Pearson_Correlation_Netatmo_DWD'] > 0.0]
#
# in_df2_after = in_df2_after[
#     in_df2_after['Bool_Pearson_Correlation_Netatmo_DWD'] > 0.0]

# in_df1 = in_df1[in_df1['Bool_Pearson_Correlation_Netatmo_DWD'] > y1_dwd0.min()]
# in_df0 = in_df0[in_df0['Bool_Pearson_Correlation_DWD_DWD'] > 0.2]
# in_df1 = in_df1[in_df1['Bool_Pearson_Correlation_DWD_DWD'] > y1_dwd0.min()]

# in_df0_dwd = in_df0_dwd[in_df0_dwd['Bool_Spearman_Correlation'] > 0.2]
# in_df1_dwd = in_df1_dwd[in_df1_dwd['Bool_Spearman_Correlation'] > 0.2]

# stns_keep_all_final_new = in_df0.index.intersection(in_df1.index)

stns_keep_all_final_new = in_df0_after.index

#cmn_stns = in_df0.index.intersection(in_df1.index)
x0_after = in_df0_after.loc[
    stns_keep_all_final_new, 'Distance to neighbor'].values.ravel()

#
# x1_after = in_df1_after.loc[
#     :, 'Distance to neighbor'].values.ravel()
#
# x2_after = in_df2_after.loc[
#     :, 'Distance to neighbor'].values.ravel()

y0_after = in_df0_after.loc[stns_keep_all_final_new,
                            'Bool_Pearson_Correlation_Netatmo_DWD'].values.ravel()

#
# y1_after = in_df1_after.loc[:,
#                             'Bool_Pearson_Correlation_Netatmo_DWD'].values.ravel()
#
#
# y2_after = in_df2_after.loc[:,
#                             'Bool_Pearson_Correlation_Netatmo_DWD'].values.ravel()

y0_after_dwd = in_df0_after.loc[stns_keep_all_final_new,
                                'Bool_Pearson_Correlation_DWD_DWD'].values.ravel()
# y1 = in_df1.loc[stns_keep_all_final_new,
#                 'Bool_Pearson_Correlation_Netatmo_DWD'].values.ravel()

# y0_dwd = in_df0.loc[:,
#                     'Bool_Pearson_Correlation_DWD_DWD'].values.ravel()
# y1_dwd = in_df1.loc[:,
#                     'Bool_Pearson_Correlation_DWD_DWD'].values.ravel()


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


stns_keep_al_sr = pd.DataFrame(data=stns_keep_all_final_new,
                               columns=['Stations'])
#
stns_keep_al_sr.to_csv(
    (save_dir /
        (r'keep_stns_all_neighbor_%s_per_%s_s0_1st_rh.csv'
         % (percent, time_freq))),
    sep=';')


# s0, x0, y0, in_df0 = read_filter_df_corr_return_stns_x_y_vals(df0)
# =============================================================================
# max_x = 34
max_x = max(x0_b4.max(), 2e4)  # x1_b4.max(), x2_b4.max()
#
plt.ioff()

fig, axs = plt.subplots(1, 2, sharex=True, sharey=True,
                        figsize=(24, 10), dpi=300)

# plt.figure(figsize=(12, 8), dpi=300)

axs[0].scatter(x0_b4, y0_b4, c='r', alpha=0.75, marker='x', s=34)
# axs[0].scatter(x1_b4, y1_b4, c='b', alpha=0.75, marker='x', s=34)
# axs[0].scatter(x2_b4, y2_b4, c='g', alpha=0.75, marker='x', s=34)

#label='First Neighbor %d Pairs' % y0.shape[0],

axs[1].scatter(x0_after, y0_after, c='r', alpha=0.75, marker='x', s=34)
# axs[1].scatter(x1_after, y1_after, c='b', alpha=0.75, marker='x', s=34)
# axs[1].scatter(x2_after, y2_after, c='g', alpha=0.75, marker='x', s=34)

axs[1].scatter(x0_dwd0, y0_dwd0, c='k', alpha=0.5, marker='o', s=34)
axs[1].scatter(x1_dwd0, y1_dwd0, c='k', alpha=0.5, marker='o', s=34)
# axs[1].scatter(x2_dwd0, y2_dwd0, c='k', alpha=0.5, marker='o', s=34)

# axs[1].scatter(x3_dwd0, y3_dwd0, c='k', alpha=0.5, marker='o', s=34)
# axs[1].scatter(x4_dwd0, y4_dwd0, c='k', alpha=0.5, marker='o', s=34)
# axs[1].scatter(x5_dwd0, y5_dwd0, c='k', alpha=0.5, marker='o', s=34)
#            label='DWD First Neighbor',
# plt.scatter(x1, y1, c='blue', alpha=0.5, marker='.',
#             label='Second Neighbor %d Pairs' % y1.shape[0], s=34)


# plt.scatter(x0, y0_dwd, c='orange', alpha=0.35, marker='.',
#             label='First Neighbor %d Pairs' % y0.shape[0], s=24)
#
# # plt.scatter(x0, y0_dwd, c='orange', alpha=0.5, marker='o',
# #             label='DWD First Neighbor', s=34)
#

# plt.scatter(x2, y2, c='green', alpha=0.35, marker='d',
#             label='Third Neighbor Stn nbr %d' % y2.shape[0], s=34)

# plt.scatter(x1, y1_dwd, c='g', alpha=0.5, marker='*',
#             label=' DWD Second Neighbor', s=34)

# plt.scatter(x0_dwd0, y0_dwd0, c='k', alpha=0.75, marker='.',
#             label='DWD First Neighbor', s=34)
# #
# plt.scatter(x1_dwd0, y1_dwd0, c='k', alpha=0.75, marker='.',
#             label='DWD Second Neighbor', s=34)


axs[0].legend(title='%d pairs' % y0_b4.shape[0], loc='upper right',  # a) left
              frameon=False, fontsize=26)._legend_box.align = 'right'
axs[1].legend(title='%d pairs' % y0_after.shape[0], loc='upper right',
              frameon=False, fontsize=26)._legend_box.align = 'right'
# ,
#             label=' DWD Third Neighbor'
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
axs[0].set_xlim([0, max_x + 300])
axs[1].set_xlim([0, max_x + 300])

axs[0].set_xticks(np.arange(0, max_x + 300, 5000))
axs[1].set_xticks(np.arange(0, max_x + 300, 5000))

# plt.xticks(np.arange(0, 21000 + 100, 5000))  # x1.max()
axs[0].set_ylim([-0.1, 1.1])
axs[0].set_xlabel('Distance [m]', labelpad=14)
axs[1].set_xlabel('Distance [m]', labelpad=14)
axs[0].set_ylabel('Indicator Correlation', labelpad=16)
# plt.legend(loc=0)
axs[0].grid(alpha=.25)
axs[1].grid(alpha=.25)
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
            (r'_%s_%s_%s_per_indic_corr_freq_%s_%s_.png'
             % (data_source0, data_source, percent, time_freq, save_acc)),
            papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
# plt.show()
