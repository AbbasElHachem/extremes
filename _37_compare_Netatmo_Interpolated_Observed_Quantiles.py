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

path_to_data = main_dir / r'oridinary_kriging_compare_DWD_Netatmo/needed_dfs'

# DAILY DATA
path_to_netatmo_interpolated_warm_season_cdf_daily_data = (
    path_to_data /
    r'interpolated_daily_data_from_qunatiles_warm_season.csv')
path_to_netatmo_interpolated_cold_season_cdf_daily_data = (
    path_to_data /
    r'interpolated_daily_data_from_qunatiles_cold_season.csv')


path_to_dwd_interpolated_cold_season_cdf_daily_data = (
    path_to_data /
    r'interpolated_dwd_daily_data_from_qunatiles_cold_season.csv')

# Original Daily Data
path_to_netatmo_daily_data = path_to_data / r'all_netatmo_ppt_data_daily_.csv'


path_to_netatmo_coords = path_to_data / r'netatmo_bw_1hour_coords_utm32.csv'

# NETATMO FIRST FILTER
path_to_netatmo_gd_stns = (path_to_data /
                           r'keep_stns_all_neighbor_90_per_60min_.csv')

# DAILY DATA
path_to_dwd_daily_data = (path_to_data /
                          r'all_dwd_daily_ppt_data_combined_2014_2019_.csv')


# HOURLY DATA
path_to_dwd_hourly_data = (path_to_data /
                           r'all_dwd_hourly_ppt_data_combined_2014_2019_.csv')

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


if use_hourly_data:
    path_to_dwd_ppt_data = path_to_dwd_hourly_data

#==============================================================================
#
#==============================================================================
netatmo_in_coords_df = pd.read_csv(path_to_netatmo_coords,
                                   index_col=0,
                                   sep=';',
                                   encoding='utf-8')

netatmo_in_vals_df = pd.read_csv(
    path_to_netatmo_daily_data, sep=';',
    index_col=0,
    encoding='utf-8', engine='c')

netatmo_in_vals_df.index = pd.to_datetime(
    netatmo_in_vals_df.index, format=idx_time_fmt)

netatmo_in_vals_df = netatmo_in_vals_df.loc[strt_date:end_date, :]
netatmo_in_vals_df.dropna(how='all', axis=0, inplace=True)
# daily sums
netatmo_in_vals_df = netatmo_in_vals_df[(0 <= netatmo_in_vals_df) &
                                        (netatmo_in_vals_df <= 300)]

cmn_stns = netatmo_in_coords_df.index.intersection(netatmo_in_vals_df.columns)
netatmo_in_vals_df = netatmo_in_vals_df.loc[:, cmn_stns]
#==============================================================================

if use_netatmo_gd_stns:

    df_gd_stns = pd.read_csv(path_to_netatmo_gd_stns,
                             index_col=0,
                             sep=';',
                             encoding='utf-8')
    good_netatmo_stns = df_gd_stns.loc[:, 'Stations'].values.ravel()
    in_vals_df = netatmo_in_vals_df.loc[:, good_netatmo_stns]
    netatmo_in_coords_df = netatmo_in_coords_df.loc[good_netatmo_stns, :]
    cmn_stns = netatmo_in_coords_df.index.intersection(
        netatmo_in_vals_df.columns)
    netatmo_in_vals_df = netatmo_in_vals_df.loc[:, cmn_stns]
#==============================================================================
# DWD DATA
#==============================================================================

dwd_in_vals_df = pd.read_csv(
    path_to_dwd_ppt_data, sep=';', index_col=0, encoding='utf-8')
dwd_in_vals_df.index = pd.to_datetime(
    dwd_in_vals_df.index, format=idx_time_fmt)

dwd_in_vals_df = dwd_in_vals_df.loc[strt_date:end_date, :]
dwd_in_vals_df.dropna(how='all', axis=0, inplace=True)


#==============================================================================
# READ DATA
#==============================================================================

df_interpolated_warm = pd.read_csv(
    path_to_netatmo_interpolated_warm_season_cdf_daily_data,
    sep=';', index_col=0)

df_interpolated_cold = pd.read_csv(
    path_to_netatmo_interpolated_cold_season_cdf_daily_data,
    sep=';', index_col=0)

df_dwd_interpolated_cold = pd.read_csv(
    path_to_dwd_interpolated_cold_season_cdf_daily_data,
    sep=';', index_col=0)

# select Netamo warm and cold data
netatmo_ppt_data_warm_season = select_season(netatmo_in_vals_df,
                                             warm_season_month)
netatmo_ppt_data_cold_season = select_season(netatmo_in_vals_df,
                                             cold_season_month)


#==============================================================================
# COMPARE
#==============================================================================
if do_it_warm_season:
    interpolated_df_season = df_interpolated_warm
    netatmo_ppt_data_season = netatmo_ppt_data_warm_season
    title_add = 'Mai_till_September_Warm'
    df_dwd_distriubutions_season = select_season(
        dwd_in_vals_df, warm_season_month)

if do_it_cold_season:
    interpolated_df_season = df_interpolated_cold

    interpolated_dwd_df_season = df_dwd_interpolated_cold
    netatmo_ppt_data_season = netatmo_ppt_data_cold_season
    title_add = 'October_till_April_Cold'

    df_dwd_distriubutions_season = select_season(
        dwd_in_vals_df, cold_season_month)


for _netatmo_stn in interpolated_df_season.columns:
    print('station is ', _netatmo_stn)

    try:
        interpolated_vals = interpolated_df_season.loc[:, _netatmo_stn]

        ppt_vals_observed = netatmo_ppt_data_season.loc[:, _netatmo_stn].dropna(
            how='all')

        ppt_vals_season, edf_vals_season = build_edf_fr_vals(
            ppt_vals_observed.values)

        plt.ioff()
        fig = plt.figure(figsize=(16, 12), dpi=150)

        ax = fig.add_subplot(111)

        # plot the stations in shapefile, look at the results of agreements

        ax.scatter(ppt_vals_season,
                   edf_vals_season,
                   alpha=.8,
                   c='b',  # colors_arr,
                   s=15,
                   marker='d',
                   # cmap=plt.get_cmap('viridis'),
                   label='Obseved Netatmo values')

        ax.scatter(interpolated_vals.values[:-1],
                   interpolated_vals.index[:-1],
                   alpha=.95,
                   c='r',  # colors_arr,
                   s=15,
                   marker='d',
                   # cmap=plt.get_cmap('viridis'),
                   label='Interpolated Netatmo values')

        for dwd_stn in df_dwd_distriubutions_season.columns:
            # print('DWD STN is', dwd_stn)
            dwd_df_stn = df_dwd_distriubutions_season.loc[:, dwd_stn].dropna()
            ppt_dwd_season, edf_dwd_season = build_edf_fr_vals(
                dwd_df_stn.values)

            ax.scatter(ppt_dwd_season,
                       edf_dwd_season,
                       alpha=.025,
                       c='k',  # colors_arr,
                       s=4,
                       marker='.')

#         dwd_interpolated_vals = interpolated_dwd_df_season.loc[:, '00020']

#         ax.scatter(dwd_interpolated_vals.values,
#                    dwd_interpolated_vals.index,
#                    alpha=.85,
#                    c='g',  # colors_arr,
#                    s=10,
#                    marker='.',
#                    label='interpolated DWD')

        ax.set_title('%s season Netatmo station is %s'
                     % (title_add, _netatmo_stn))
        ax.grid(alpha=0.25)
        ax.set_ylim([.5, 1.01])
        ax.set_xlim([-0.1, max(interpolated_vals.values.max(),
                               ppt_vals_season.max()) + 2])
        ax.set_xlabel('Rainfall (mm/day)')
        ax.legend(loc='lower right')
        ax.set_ylabel('CDF')

        plt.savefig(out_plots_path / (
            'oberserved_vs_interpolated_cdf_%s_season_%s_.png'
            % (title_add, _netatmo_stn)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
        plt.close()

    except Exception as msg:
        print(msg)
        continue
