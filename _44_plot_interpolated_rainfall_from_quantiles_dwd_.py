'''
Created on 4 Sep 2019

@author: hachem
'''

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from pathlib import Path
from _00_additional_functions import (build_edf_fr_vals, select_season)

plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})

# =============================================================================

main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')
os.chdir(main_dir)

out_plots_path = main_dir / r'oridinary_kriging_compare_DWD_Netatmo'

path_to_data = main_dir / r'NetAtmo_BW'

# DAILY DATA
path_to_dwd_daily_data = (path_to_data /
                          r'all_dwd_daily_ppt_data_combined_2014_2019_.csv')

# # HOURLY DATA
# path_to_dwd_hourly_data = (path_to_data /
#                            r'all_dwd_hourly_ppt_data_combined_2014_2019_.csv')

#==============================================================================
#
#==============================================================================
strt_date = '2014-01-01'
end_date = '2019-08-01'
idx_time_fmt = '%Y-%m-%d'


warm_season_month = [5, 6, 7, 8, 9]  # mai till sep
cold_season_month = [10, 11, 12, 1, 2, 3, 4]  # oct till april

use_netatmo_gd_stns = True

use_daily_data = True  # False
use_hourly_data = False  # True

do_it_warm_season = True  # False
do_it_cold_season = False  # True

if use_daily_data:
    path_to_dwd_ppt_data = path_to_dwd_daily_data

# if use_hourly_data:
#     path_to_dwd_ppt_data = path_to_dwd_hourly_data
#==============================================================================
# DWD DATA
#==============================================================================

dwd_in_vals_df = pd.read_csv(
    path_to_dwd_ppt_data, sep=';', index_col=0, encoding='utf-8')
dwd_in_vals_df.index = pd.to_datetime(
    dwd_in_vals_df.index, format=idx_time_fmt)

dwd_in_vals_df = dwd_in_vals_df.loc[strt_date:end_date, :]
dwd_in_vals_df.dropna(how='all', axis=0, inplace=True)

season = 'warm'
title_add = 'Mai_till_September_Warm'  # title_add = 'October_till_April_Cold'

df_dwd_distriubutions_season = select_season(
    dwd_in_vals_df, warm_season_month)

for i in range(1):

    path_to_dwd_quantile_interpolated_using_dwd_netatmo_season_cdf_daily_data = (
        r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
        r"\oridinary_kriging_compare_DWD_Netatmo"
        r"\interpolated_ppt_dwd_daily_data_basedon_quantiles_%s_season_using_dwd_netamo.csv"
        % (season))

    path_to_dwd_quantile_interpolated_using_dwd_season_cdf_daily_data = (
        r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
        r"\oridinary_kriging_compare_DWD_Netatmo"

        r"\interpolated_ppt_dwd_daily_data_basedon_qunatiles_%s_season_using_dwd_only.csv"
        % (season))

    path_to_dwd_quantile_interpolated_using_netatmo_season_cdf_daily_data = (
        r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
        r"\oridinary_kriging_compare_DWD_Netatmo"
        r"\interpolated_dwd_daily_data_basedon_qunatiles_%s_season_using_netatmo_only.csv"
        % (season))

    #=========================================================================
    # INTERPOLATION OF QUANTILES BASED ON RAINFALL
    #=========================================================================
    # COLD SEASON
    # DWD stations interpolated using DWD and Netatmo stns
    df_dwd_ppt_interpolated_using_dwd_netatmo = pd.read_csv(
        path_to_dwd_quantile_interpolated_using_dwd_netatmo_season_cdf_daily_data,
        sep=';', index_col=0)

    # DWD stations interpolated using DWD stns
    df_dwd_ppt_interpolated_using_dwd_only = pd.read_csv(
        path_to_dwd_quantile_interpolated_using_dwd_season_cdf_daily_data,
        sep=';', index_col=0)

    # DWD stations interpolated using Netatmo stns
    df_dwd_ppt_interpolated_using_netatmo_only = pd.read_csv(
        path_to_dwd_quantile_interpolated_using_netatmo_season_cdf_daily_data,
        sep=';', index_col=0)

    for _dwd_stn_ in df_dwd_ppt_interpolated_using_dwd_netatmo.columns:
        print('station is ', _dwd_stn_)

        try:

            # interpolated quantiles based on quantiles
            int_ppt_vals_dwd_netatmo_used = df_dwd_ppt_interpolated_using_dwd_netatmo.loc[
                :, _dwd_stn_]
            int_ppt_vals_dwd_used = df_dwd_ppt_interpolated_using_dwd_only.loc[
                :, _dwd_stn_]
            int_ppt_vals_dwd_netatmo_only = df_dwd_ppt_interpolated_using_netatmo_only.loc[
                :, _dwd_stn_]

            # oroginal data
            ppt_vals_observed = df_dwd_distriubutions_season.loc[:, _dwd_stn_].dropna(
                how='all')
            ppt_vals_season, edf_vals_season = build_edf_fr_vals(
                ppt_vals_observed.values)

            plt.ioff()
            fig = plt.figure(figsize=(16, 12), dpi=150)

            ax = fig.add_subplot(111)

            # plot the stations in shapefile, look at the results of agreements

            ax.scatter(ppt_vals_season,
                       edf_vals_season,
                       alpha=.2,
                       c='grey',  # colors_arr,
                       s=10,
                       marker='d',
                       # cmap=plt.get_cmap('viridis'),
                       label='Obseved DWD values')

            #==================================================================
            # PLOT INTERPOLATED QUANTILES BASED ON PPT
            #======================================================================,

            ax.scatter(
                int_ppt_vals_dwd_used.index,
                int_ppt_vals_dwd_used.values,
                alpha=.75,
                c='r',  # colors_arr,
                s=15,
                marker='o',
                # cmap=plt.get_cmap('viridis'),
                       label='Interpolated DWD Quantiles values with DWD stns')
            ax.scatter(
                int_ppt_vals_dwd_netatmo_used.index,
                int_ppt_vals_dwd_netatmo_used.values,
                alpha=.75,
                c='b',  # colors_arr,
                s=15,
                marker='+',
                # cmap=plt.get_cmap('viridis'),
                       label='Interpolated DWD Quantiles values with DWD & Netatmo stns')

            ax.scatter(
                int_ppt_vals_dwd_netatmo_only.index,
                int_ppt_vals_dwd_netatmo_only.values,
                alpha=.75,
                c='g',  # colors_arr,
                s=15,
                marker='X',
                # cmap=plt.get_cmap('viridis'),
                       label='Interpolated DWD Quantiles values with Netatmo stns')

            ax.set_title('%s season DWD station is %s'
                         % (title_add, _dwd_stn_))
            ax.grid(alpha=0.25)
            ax.set_ylim([0.5, 1.01])
    #         ax.set_xlim([-0.5, max(int_vals_dwd_netatmo_used.values.max(),
    #                                int_vals_dwd_used.max(),
    #                                int_vals_dwd_netatmo_only.values.max(),
    #                                ppt_vals_season.max()) + 2])
            ax.set_xlabel('Rainfall (mm/day)')
            ax.legend(loc='lower right')
            ax.set_ylabel('CDF')

            plt.savefig(out_plots_path / (
                'oberserved_vs_interpolated_cdf_from_ppt_%s_season_%s.png'
                % (title_add, _dwd_stn_)),
                frameon=True, papertype='a4',
                bbox_inches='tight', pad_inches=.2)
            plt.close()

        except Exception as msg:
            print(msg)
            continue
