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


from _00_additional_functions import (select_df_within_period,
                                      calculate_probab_ppt_below_thr,
                                      build_edf_fr_vals, resampleDf)

plt.ioff()
# =============================================================================

# def values to replace edf of ppt=0
ppt_min_thr_0_vals = 0.1  # everything below it gets value of P0

netatmo_data = False
use_good_netatmo_stns = False

dwd_data = True
resample_data = True

# select data only within this period (same as netatmo / dwd)
start_date = '2015-01-01 00:00:00'
end_date = '2019-09-01 00:00:00'

resample_frequencies = ['60min', '120min', '180min',
                        '360min', '720min', '1440min']
# =============================================================================

if netatmo_data:
    #     df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    #                r"\all_netatmo_ppt_data_monthly_.csv")
    #     df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    #                r"\all_netatmo_ppt_data_daily_.csv")
    df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
               # r"\all_netatmo_ppt_data_daily_.csv")
               r"\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.csv")
    save_acc = 'netatmo'
# list of Netatmo stations with good indicator correlations
if use_good_netatmo_stns:
    df_stns = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
               r"\plots_NetAtmo_ppt_DWD_ppt_correlation_"
               r"\keep_stns_all_neighbor_99_0_per_60min_s0.csv")
    save_acc = 'netatmo_good_stns'
if dwd_data:
    df_file = (
        # \all_dwd_ppt_data_daily_.csv")
        r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
        r"\all_dwd_hourly_ppt_data_combined_2014_2019_.csv")

    df_file = (
        # \all_dwd_ppt_data_daily_.csv")
        r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
        r"\ppt_all_dwd_old_60min_.csv")

    save_acc = 'dwd'
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

# in_df = select_df_within_period(in_df, start=start_date, end=end_date)

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


if resample_data:

    for agg_freq in resample_frequencies:

        print('Stations resampling for', agg_freq)
        df_stn0 = resampleDf(df=in_df,
                             agg=agg_freq)

        all_stns = df_stn0.columns.to_list()

        data_mtx = np.zeros(shape=(df_stn0.index.shape[0],
                                   len(all_stns))).astype('float')
        data_mtx[data_mtx == 0] = np.nan

        df_all = pd.DataFrame(index=df_stn0.index,
                              columns=all_stns,
                              data=data_mtx)
        for stn_ in df_stn0.columns:
            #print('stn is ', stn_)
            # %% transform data to quantiles from cdf
            try:
                stn_df_no_nans = df_stn0.loc[:, stn_].dropna()
                if stn_df_no_nans.size > 2:
                    p0 = calculate_probab_ppt_below_thr(stn_df_no_nans.values,
                                                        ppt_min_thr_0_vals)
                    # print('P0', p0, p0 / 2)
                    x0, y0 = build_edf_fr_vals(stn_df_no_nans.values)

                    y0[np.where(x0 <= 0.1)] = p0 / 2
                    # TODO: Check again
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

        df_all.dropna(how='all', inplace=True)
        df_all.to_csv((r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
                       r"\edf_ppt_all_%s_old_%s_.csv" % (save_acc, agg_freq)),
                      sep=';', float_format='%.3f')

    print('DONE WITH EVERYTHNG !')


raise Exception


for stn_ in df_stn0.columns:
    print('stn is ', stn_)
    # %% transform data to quantiles from cdf
    try:
        stn_df_no_nans = df_stn0.loc[:, stn_].dropna()
        if stn_df_no_nans.size > 2:
            p0 = calculate_probab_ppt_below_thr(stn_df_no_nans.values,
                                                ppt_min_thr_0_vals)
            print('P0', p0, p0 / 2)
            x0, y0 = build_edf_fr_vals(stn_df_no_nans.values)

            y0[np.where(x0 <= 0.1)] = p0 / 2

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
               r"\edf_ppt_all_%s_hourly_.csv" % save_acc),
              sep=';', float_format='%.3f')

print('DONE WITH EVERYTHNG !')
