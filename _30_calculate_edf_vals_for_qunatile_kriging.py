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

from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from _00_additional_functions import build_edf_fr_vals

plt.ioff()
#==============================================================================

# def values to replace edf of ppt=0
# TODO: TALK TO PROF
edf_ppt_0 = 0.44

# create indicator Kriging

df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
           r"\ppt_all_netatmo_hourly_stns_combined_.fk")

# list of stations with good indicator correlations

df_stns = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
           r"\filter_Netamo_data_basedon_indicator_correlation"
           r"\keep_stns_all_neighbors_combined_95_per_.csv")

in_df = pd.read_feather(df_file, use_threads=True)
in_df.set_index('Time', inplace=True)
in_df.index = pd.to_datetime(in_df.index, format='%Y-%m-%d %H:%M:%S')
in_df.dropna(how='all', inplace=True)
# all_stns = in_df.columns.to_list()
#==============================================================================

in_df_stns = pd.read_csv(df_stns, index_col=0, sep=';', header=None)
good_stns = in_df_stns.values.ravel()

#==============================================================================
# create empty dataframe to hold the results
data_mtx = np.zeros(
    shape=(in_df.index.shape[0], len(good_stns))).astype('float')
data_mtx[data_mtx == 0] = np.nan

df_all = pd.DataFrame(index=in_df.index, columns=good_stns, data=data_mtx)

#==============================================================================
# select from the original station the staions with 'good' data
df_stn0 = in_df.loc[:, good_stns].dropna(how='all')


for stn_ in df_stn0.columns:
    print('stn is ', stn_)
    # %% transform data to quantiles from cdf
    stn_df_no_nans = df_stn0.loc[:, stn_].dropna()

    x0, y0 = build_edf_fr_vals(stn_df_no_nans.values)
    y0[np.where(x0 == 0)] = edf_ppt_0

    df_stn_ = pd.DataFrame(data=stn_df_no_nans.values,
                           index=stn_df_no_nans.index,
                           columns=[stn_])
    df_stn_.sort_values(by=stn_, inplace=True)

    df_stn_.loc[:, 'edf'] = y0
    df_stn_.sort_index(inplace=True)
    df_all.loc[:, stn_] = df_stn_.loc[:, 'edf']

df_all.dropna(how='all', inplace=True)
df_all.to_csv((r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
               r"\edf_ppt_all_netatmo_hourly_stns_combined_.csv"),
              sep=';', float_format='%.3f')

print('DONE WITH EVERYTHNG !')
# # %%
#
# fig = plt.figure(figsize=(12, 8), dpi=100)
# plt.plot(x0, y0, c='r', alpha=.5)

# plt.show()
#%%
