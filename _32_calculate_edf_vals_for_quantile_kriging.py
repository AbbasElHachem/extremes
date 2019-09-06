# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Convert PPT values to EDF values for Quantile Kriging
Purpose: Prepare data to be used for quantile kriging

Created on: 2019-08-12

Parameters
----------
    Input: Rainfall data

Returns
-------
    Output: Df with index as time and values as EDF for every station
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

# TODO: make soft coded

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from _00_additional_functions import build_edf_fr_vals

plt.ioff()
# =============================================================================

# def values to replace edf of ppt=0
# TODO: TALK TO PROF
edf_ppt_0 = 0.44

netatmo_data = True
use_good_netatmo_stns = False

dwd_data = False

# =============================================================================

if netatmo_data:
    #     df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    #                r"\all_netatmo_ppt_data_monthly_.csv")
    #     df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    #                r"\all_netatmo_ppt_data_daily_.csv")
    df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
               r"\ppt_all_netatmo_hourly_stns_combined_new.csv")
# list of Netatmo stations with good indicator correlations
if use_good_netatmo_stns:
    df_stns = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
               r"\plots_NetAtmo_ppt_DWD_ppt_correlation_"
               r"\keep_stns_all_neighbor_92_per_60min_.csv")

if dwd_data:
    df_file = (r"F:\download_DWD_data_recent"
               r"\all_dwd_hourly_ppt_data_combined_2014_2019_.csv")
    # df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    #            r"\all_dwd_ppt_data_monthly_.csv")


#==============================================================================
# # read DF ppt data
#==============================================================================
in_df = pd.read_csv(df_file, index_col=0, sep=';',
                    infer_datetime_format=True,
                    parse_dates=True, engine='c')

in_df.index = pd.to_datetime(in_df.index,
                             format='%Y-%m-%d %H:%M:%S')
in_df.dropna(how='all', inplace=True)

all_stns = in_df.columns.to_list()
#==============================================================================
# create empty dataframe to hold the results
if not use_good_netatmo_stns:
    data_mtx = np.zeros(shape=(in_df.index.shape[0],
                               len(all_stns))).astype('float')
    data_mtx[data_mtx == 0] = np.nan

    df_all = pd.DataFrame(index=in_df.index,
                          columns=all_stns,
                          data=data_mtx)
    df_stn0 = in_df.loc[:, all_stns].dropna(how='all')

#==============================================================================
# select from the original station the staions with 'good' data

if use_good_netatmo_stns:
    in_df_stns = pd.read_csv(df_stns, index_col=0, sep=';')
    good_stns = in_df_stns.values.ravel()

    data_mtx = np.zeros(
        shape=(in_df.index.shape[0], len(good_stns))).astype('float')
    data_mtx[data_mtx == 0] = np.nan

    df_all = pd.DataFrame(index=in_df.index, columns=good_stns, data=data_mtx)
    df_stn0 = in_df.loc[:, good_stns].dropna(how='all')
#==============================================================================

for stn_ in df_stn0.columns:
    print('stn is ', stn_)
    # %% transform data to quantiles from cdf
    try:
        stn_df_no_nans = df_stn0.loc[:, stn_].dropna()
        if stn_df_no_nans.size > 2:
            x0, y0 = build_edf_fr_vals(stn_df_no_nans.values)
            y0[np.where(x0 == 0)] = edf_ppt_0

            df_stn_ = pd.DataFrame(data=stn_df_no_nans.values,
                                   index=stn_df_no_nans.index,
                                   columns=[stn_])
            df_stn_.sort_values(by=stn_, inplace=True)

            df_stn_.loc[:, 'edf'] = y0
            df_stn_.sort_index(inplace=True)
            df_all.loc[:, stn_] = df_stn_.loc[:, 'edf']
        else:
            df_all.loc[:, stn_] = np.nan
    except Exception as msg:
        print(msg)
        continue
        #     break

df_all.dropna(how='all', inplace=True)
df_all.to_csv((r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
               r"\edf_ppt_all_netatmo_hourly_.csv"),
              sep=';', float_format='%.3f')

print('DONE WITH EVERYTHNG !')
