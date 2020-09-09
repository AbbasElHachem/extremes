# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Calculate and plot statistical differences between neighbours
Purpose: Find validity of Netatmo Station compared to DWD station

"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

import os
import timeit
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib import rc
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears
from scipy.spatial import cKDTree

from _00_additional_functions import (resample_intersect_2_dfs,
                                      select_convective_season,
                                      select_df_within_period,
                                      build_edf_fr_vals,
                                      get_cdf_part_abv_thr,
                                      plt_on_map_agreements,
                                      plt_correlation_with_distance)


#==============================================================================
#
#==============================================================================
main_dir = Path(os.getcwd())
os.chdir(main_dir)

register_matplotlib_converters()

plt.ioff()

rc('font', size=16)
rc('font', family='arial')
rc('axes', labelsize=16)
rcParams['axes.labelpad'] = 10
plt.rcParams['axes.labelweight'] = 'bold'
#==============================================================================
# # for getting station names in BW
#==============================================================================
path_to_ppt_netatmo_data_csv = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.csv')
assert os.path.exists(path_to_ppt_netatmo_data_csv), 'wrong NETATMO Ppt file'

#==============================================================================
# # for reading ppt data station by station
# # HOURLY DATA
#==============================================================================

path_to_ppt_netatmo_data_feather = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.fk')

assert os.path.exists(
    path_to_ppt_netatmo_data_feather), 'wrong NETATMO Ppt file'

# HOURLY DATA

path_to_ppt_dwd_data = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\all_dwd_hourly_ppt_data_combined_2015_2019_.fk")
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'

# BW

path_to_netatmo_coords_df_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\rain_bw_1hour"
    r"\netatmo_bw_1hour_coords.csv")
assert os.path.exists(path_to_netatmo_coords_df_file), 'wrong DWD coords file'

#==============================================================================
# # coords of stns in utm32 in BW
#==============================================================================

# COORDINATES
path_to_dwd_coords_utm32 = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r'/station_coordinates_names_hourly_only_in_BW_utm32.csv')

path_to_netatmo_coords_utm32 = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r'/netatmo_bw_1hour_coords_utm32.csv')

out_save_dir_orig = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                     r'\indic_corr_2019')

if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

# path to remaining stns
path_to_remaining_stns = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\indic_corr_2019"
    r"\keep_stns_all_neighbor_99_0_per_60min_s0_1st_rh.csv"
)

df0_b4 = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\indic_corr_2019"
    r"\2pearson_year_allyears_df_comparing_correlations_max_sep_dist_2500_freq_60min_dwd_netatmo_upper_99_0_percent_data_considered_neighbor_0_.csv"
)

in_df0_b4 = pd.read_csv(df0_b4, index_col=0, sep=';').dropna(how='all')

#==============================================================================
#
#==============================================================================
# used in Netatmo coords df
x_col_name = 'lon'
y_col_name = 'lat'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 2500  # 5000  # m

# threshold for max ppt value per hour
max_ppt_thr = 200.  # ppt above this value are not considered

# only highest x% of the values are selected
lower_percentile_val_lst = [99.0]  # [80, 85, 90, 95, 99]

# ['10min', '60min', '120min', '480min', '720min', '1440min']
aggregation_frequencies = ['60min']

# temporal aggregation of df

# [0, 1, 2, 3, 4]  # refers to DWD neighbot (0=first)
neighbors_to_chose_lst = [0]
# 30 days * 24 hours * 2month
# minimum hourly values that should be available per station
min_req_ppt_vals = 10  # 30 * 24 * 2

# this is used to keep only data where month is not in this list
# not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april
not_convective_season = []  # oct till april

date_fmt = '%Y-%m-%d %H:%M:%S'

# select data only within this period (same as netatmo)
start_date = '2019-01-01 00:00:00'
end_date = '2019-12-30 00:00:00'

nbr_years = 1  # 2015 +4 = 2019
years = ['2019']  # , '2016', '2017', '2018', '2019']
#==============================================================================
#
#==============================================================================


def plot_histogram(x, density=True, nbins=10):
    hist, bins = np.histogram(x, bins=nbins, density=density)
    width = 0.5 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    return hist, width, center


# after
netatmo_df_after = pd.read_csv(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\indic_corr_2019\new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_2500_freq_60min_dwd_netatmo_upper_99_0_percent_data_considered_neighbor_0_.csv",
    sep=';', index_col=0)
# netatmo ids remaining
netatmo_df_after_filter = pd.read_csv(
    path_to_remaining_stns, index_col=1, sep=';')
netatmo_ids_after_filter = netatmo_df_after_filter.index.to_list()

neatmo_ids_b4_filter = in_df0_b4.index.to_list()
not_filtered_ids = list(set(neatmo_ids_b4_filter) -
                        set(netatmo_ids_after_filter))

in_netatmo_df_coords_utm32 = pd.read_csv(
    path_to_netatmo_coords_utm32, sep=';',
    index_col=0, engine='c')
netatmo_df_after.columns, in_df0_b4.columns
ranks_b4 = in_df0_b4.loc[not_filtered_ids, 'Orig_Spearman_Correlation']
ranks_b4 = ranks_b4[ranks_b4 < 1]
i1 = np.where(
    netatmo_df_after.loc[:, 'Bool_Pearson_Correlation_Netatmo_DWD'].values < 1)[0]
i0 = np.where(
    netatmo_df_after.loc[:, 'Bool_Pearson_Correlation_Netatmo_DWD'].values > 0.2)[0]
i10 = np.intersect1d(i1, i0)
ix_ = netatmo_df_after.index[i10]
ranks_after = netatmo_df_after.loc[ix_,
                                   'Orig_Spearman_Correlation_Netatmo_DWD']
# ranks_after = ranks_after[ranks_after < 1]
# ranks_after = ranks_after[ranks_after > 0.2]
ranks_after.min()
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
plt.ioff()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 12),
                               sharex=True, sharey=True, dpi=100)
hist0, width0, center0 = plot_histogram(
    ranks_b4.values.ravel(), density=False, nbins=10)


hist1, width1, center1 = plot_histogram(
    ranks_after.values.ravel(), density=False, nbins=10)
max_y = max(max(hist0), max(hist1))
ax1.set_ylim([0, max_y + 3])
ax2.set_yticks(np.arange(0, max_y + 3, 5))

ax1.bar(center0, hist0,  align='center', width=0.068,  # 0.0448
        alpha=0.95, color='lime', edgecolor='g',
        linewidth=1, label='n=%d' % len(ranks_b4))

ax1.set_xlabel('Rank Correlation', fontsize=18)
ax1.grid(alpha=0.35)
ax1.set_ylabel('Frequency', fontsize=18)
ax1.legend(loc='upper right')
ax1.title.set_text('Removed PWS ')
ax1.set_xlim([0, 1])
ax1.set_xticks(np.arange(0, 1.01, .1))


ax1.xaxis.set_minor_locator(AutoMinorLocator())

ax1.tick_params(which='both', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=3, color='grey')


ax2.bar(center1, hist1, align='center', width=0.0715,  # 0.0538
        alpha=0.95, color='lime', edgecolor='g',
        linewidth=1, label=' n=%d' % len(
            ranks_after))
ax2.set_xlim([0, 1])
ax2.set_xlabel('Rank Correlation', fontsize=18)  # Indicator
ax2.grid(alpha=0.35)
ax2.set_xticks(np.arange(0, 1.01, .1))
# ax2.set_ylabel('Frequency')
ax2.legend(loc='upper right')
ax2.title.set_text('Remaining PWS')

ax2.xaxis.set_minor_locator(AutoMinorLocator())

ax2.tick_params(which='both', width=2)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=3, color='grey')

fig.tight_layout()
fig.savefig(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\indic_corr_2019\corr_.png",
    frameon=True, papertype='a4',
    bbox_inches='tight', pad_inches=.2)
plt.close()
# fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(20, 12),
#                                 sharex=True, sharey=True)
#
# for ppt_stn_id in netatmo_ids_after_filter:
#     # iterating through netatmo ppt stations
#
#     print('\n********\n First Ppt Stn Id is', ppt_stn_id)
#
#     try:
#         netatmo_ppt_stn1_orig = pd.read_feather(path_to_ppt_netatmo_data_feather,
#                                                 columns=[
#                                                     'Time', ppt_stn_id],
#                                                 use_threads=True)
#
#         netatmo_ppt_stn1_orig.set_index('Time', inplace=True)
#     except Exception as msg:
#             # print('error reading dwd', msg)
#         netatmo_ppt_stn1_orig = pd.read_feather(path_to_ppt_netatmo_data_feather,
#                                                 columns=[
#                                                     'index', ppt_stn_id],
#                                                 use_threads=True)
#         netatmo_ppt_stn1_orig.set_index('index', inplace=True)
#
#     netatmo_ppt_stn1_orig.index = pd.to_datetime(
#         netatmo_ppt_stn1_orig.index, format=date_fmt)
#
#     # drop all index with nan values
#     netatmo_ppt_stn1_orig.dropna(axis=0, inplace=True)
#
#     netatmo_ppt_stn1_orig = netatmo_ppt_stn1_orig[
#         netatmo_ppt_stn1_orig < max_ppt_thr]
#
#     # select only convective season
#     netatmo_ppt_stn1_orig = select_convective_season(
#         df=netatmo_ppt_stn1_orig,
#         month_lst=not_convective_season)
#     netatmo_ppt_stn1_orig = select_df_within_period(netatmo_ppt_stn1_orig,
#                                                     start_date, end_date)
#
#     x1, y1 = build_edf_fr_vals(netatmo_ppt_stn1_orig.values.ravel())
#     ax2.plot(x1, y1, c='k', alpha=0.25)
#     ax2.set_xlabel('[mm/h]')
#     ax2.grid(alpha=0.25)
# #     ax2.legend(loc='upper right')
#     ax2.title.set_text('Remaining stations n=%d' % len(
#         netatmo_ids_after_filter))
#
#     netatmo_ppt_stn1_orig = netatmo_ppt_stn1_orig[netatmo_ppt_stn1_orig > 0]
#     netatmo_ppt_stn1_orig[netatmo_ppt_stn1_orig > 40] = 40
#     netatmo_ppt_stn1_orig.dropna(inplace=True)
#     hist1, width1, center1 = plot_histogram(
#         netatmo_ppt_stn1_orig.values.ravel(), density=False)
#
#     ax4.bar(center1, hist1, align='center', width=width1,
#             alpha=0.25, color='b', edgecolor='darkblue',
#             linewidth=1)
#
#     ax4.set_xlabel('[mm/h]')
#     ax4.grid(alpha=0.25)
# #     ax2.legend(loc='upper right')
#     ax4.title.set_text('Remaining stations n=%d' % len(
#         netatmo_ids_after_filter))
# for ppt_stn_id in not_filtered_ids:
#     # iterating through netatmo ppt stations
#
#     print('\n********\n First Ppt Stn Id is', ppt_stn_id)
#
#     try:
#         netatmo_ppt_stn1_orig = pd.read_feather(path_to_ppt_netatmo_data_feather,
#                                                 columns=[
#                                                     'Time', ppt_stn_id],
#                                                 use_threads=True)
#
#         netatmo_ppt_stn1_orig.set_index('Time', inplace=True)
#     except Exception as msg:
#             # print('error reading dwd', msg)
#         netatmo_ppt_stn1_orig = pd.read_feather(path_to_ppt_netatmo_data_feather,
#                                                 columns=[
#                                                     'index', ppt_stn_id],
#                                                 use_threads=True)
#         netatmo_ppt_stn1_orig.set_index('index', inplace=True)
#
#     netatmo_ppt_stn1_orig.index = pd.to_datetime(
#         netatmo_ppt_stn1_orig.index, format=date_fmt)
#
#     # drop all index with nan values
#     netatmo_ppt_stn1_orig.dropna(axis=0, inplace=True)
#
#     netatmo_ppt_stn1_orig = netatmo_ppt_stn1_orig[
#         netatmo_ppt_stn1_orig < max_ppt_thr]
#
#     # select only convective season
#     netatmo_ppt_stn1_orig = select_convective_season(
#         df=netatmo_ppt_stn1_orig,
#         month_lst=not_convective_season)
#     netatmo_ppt_stn1_orig = select_df_within_period(netatmo_ppt_stn1_orig,
#                                                     start_date, end_date)
#
# #         hist1, width1, center1 = plot_histogram(
# #             netatmo_ppt_stn1_orig.values.ravel(), density=True)
# #
# #         ax.bar(center1, hist1, align='center', width=width1,
# #                alpha=0.5, color='b', edgecolor='darkblue',
# #                linewidth=1)
#
#     x1, y1 = build_edf_fr_vals(netatmo_ppt_stn1_orig.values.ravel())
#     ax1.plot(x1, y1, c='k', alpha=0.25)
#     ax1.set_xlabel('[mm/h]')
#     ax1.set_ylabel('F(x)')
#     ax1.set_ylim([0.99, 1.001])
#     ax1.grid(alpha=0.25)
# #     ax1.legend(loc='upper right')
#     ax1.title.set_text('Removed stations n=%d' % len(not_filtered_ids))
#
#     netatmo_ppt_stn1_orig = netatmo_ppt_stn1_orig[netatmo_ppt_stn1_orig > 0]
#     netatmo_ppt_stn1_orig[netatmo_ppt_stn1_orig > 40] = 40
#     netatmo_ppt_stn1_orig.dropna(inplace=True)
#     hist1, width1, center1 = plot_histogram(
#         netatmo_ppt_stn1_orig.values.ravel(), density=False)
#
#     ax3.bar(center1, hist1, align='center', width=width1,
#             alpha=0.25, color='b', edgecolor='darkblue',
#             linewidth=1)
#
#     ax3.set_xlabel('[mm/h]')
#     ax3.grid(alpha=0.25)
#     ax3.set_ylabel('F(x)')
# #     ax2.legend(loc='upper right')
#     ax3.title.set_text('Removed stations n=%d' % len(
#         not_filtered_ids))
#
# fig.savefig(
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\indic_corr_2019\cdf.png",
#     frameon=True, papertype='a4',
#     bbox_inches='tight', pad_inches=.2)
#
# fig2.savefig(
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\indic_corr_2019\hist.png",
#     frameon=True, papertype='a4',
#     bbox_inches='tight', pad_inches=.2)


# plt.show()
# plt.close()
