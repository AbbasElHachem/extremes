# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s

TODO: ADD DOCUMENTATION
"""
import os
import timeit
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from pathlib import Path
from _00_additional_functions import build_edf_fr_vals

plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})

# =============================================================================

main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')
os.chdir(main_dir)

out_plots_path = main_dir / r'oridinary_kriging_compare_DWD_Netatmo'

path_to_data = main_dir / r'oridinary_kriging_compare_DWD_Netatmo/needed_dfs'

# DAILY DATA
path_to_dwd_daily_data = (path_to_data /
                          r'all_dwd_daily_ppt_data_combined_2014_2019_.csv')


# HOURLY DATA
path_to_dwd_hourly_data = (path_to_data /
                           r'all_dwd_hourly_ppt_data_combined_2014_2019_.csv')

# =============================================================================
strt_date = '2015-01-01'
end_date = '2019-09-01'

warm_season_month = [5, 6, 7, 8, 9]  # mai till sep
cold_season_month = [10, 11, 12, 1, 2, 3, 4]  # oct till april


use_daily_data = False  # False
use_hourly_data = True  # True


# =============================================================================
if use_daily_data:
    path_to_dwd_ppt_data = path_to_dwd_daily_data
    acc = 'daily'

if use_hourly_data:
    path_to_dwd_ppt_data = path_to_dwd_hourly_data
    acc = 'hourly'

#==============================================================================
# # DWD DATA
#==============================================================================


dwd_in_vals_df = pd.read_csv(
    path_to_dwd_ppt_data, sep=';', index_col=0, encoding='utf-8')
dwd_in_vals_df.index = pd.to_datetime(
    dwd_in_vals_df.index, format='%Y-%m-%d')

dwd_in_vals_df = dwd_in_vals_df.loc[strt_date:end_date, :]
dwd_in_vals_df.dropna(how='all', axis=0, inplace=True)


#==============================================================================


def select_season(df,  # df to slice, index should be datetime
                  month_lst  # list of month for convective season
                  ):
    """
    return dataframe without the data corresponding to the winter season
    """
    df = df.copy()
    df_conv_season = df[df.index.month.isin(month_lst)]

    return df_conv_season

#==============================================================================
#
#==============================================================================


print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
start = timeit.default_timer()  # to get the runtime of the program

#==============================================================================
# # ordinary kriging
#==============================================================================

df_dwd_distriubutions_cold_season = select_season(
    dwd_in_vals_df, cold_season_month)
df_dwd_distriubutions_warm_season = select_season(
    dwd_in_vals_df, warm_season_month)


for stn_id in dwd_in_vals_df.columns:
    print('station is', stn_id)
    stn_data_df = dwd_in_vals_df.loc[:, stn_id].dropna()

    stn_data_cold_season = select_season(stn_data_df, cold_season_month)
    stn_data_warm_season = select_season(stn_data_df, warm_season_month)

    ppt_cold_season, edf_cold_season = build_edf_fr_vals(
        stn_data_cold_season.values)
    ppt_warm_season, edf_warm_season = build_edf_fr_vals(
        stn_data_warm_season.values)

    #==========================================================================
    #
    #==========================================================================
    df_stn_cold = pd.DataFrame(data=stn_data_cold_season.values,
                               index=stn_data_cold_season.index,
                               columns=[stn_id])
    df_stn_cold.sort_values(by=stn_id, inplace=True)

    df_stn_cold.loc[:, 'edf'] = edf_cold_season
    df_stn_cold.sort_index(inplace=True)

    df_dwd_distriubutions_cold_season.loc[
        stn_data_cold_season.index,
        stn_id] = df_stn_cold.loc[stn_data_cold_season.index, 'edf']
    #==========================================================================
    #
    #==========================================================================
    df_stn_warm = pd.DataFrame(data=stn_data_warm_season.values,
                               index=stn_data_warm_season.index,
                               columns=[stn_id])
    df_stn_warm.sort_values(by=stn_id, inplace=True)

    df_stn_warm.loc[:, 'edf'] = edf_warm_season
    df_stn_warm.sort_index(inplace=True)

    df_dwd_distriubutions_warm_season.loc[
        stn_data_warm_season.index,
        stn_id] = df_stn_warm.loc[stn_data_warm_season.index, 'edf']

df_dwd_distriubutions_warm_season.to_csv(
    out_plots_path / ('df_dwd_distributions_warm_season_%s.csv' % acc),
    sep=';', float_format='%.2f')

df_dwd_distriubutions_cold_season.to_csv(
    out_plots_path / ('df_dwd_distributions_cold_season_%s.csv' % acc),
    sep=';', float_format='%.2f')

print('done with everything')
