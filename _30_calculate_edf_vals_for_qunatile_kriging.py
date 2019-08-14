#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:58:27 2019

@author: abbas
"""
#%%
import pandas as pd
import numpy as np

from _00_additional_functions import build_edf_fr_vals
# create indicator Kriging

df_file = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.fk"
#%%
in_df = pd.read_feather(df_file, use_threads=True)
in_df.set_index('Time', inplace=True)
in_df.index = pd.to_datetime(in_df.index, format='%Y-%m-%d %H:%M:%S')
in_df.dropna(how='all', inplace=True)
#%%

all_stns = in_df.columns.to_list()

# list of stations with good indicator correlations

df_stns = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\filter_Netamo_data_basedon_indicator_correlation\keep_stns_all_neighbors_combined_95_per_.csv"

in_df_stns = pd.read_csv(df_stns, index_col=0, sep=';', header=None)

good_stns = in_df_stns.values.ravel()

# date_range = pd.date_range('2012-01-01 00:00:00',
#                            '2019-07-31 00:00:00',
#                            freq='H')  # 'H'

data_mtx = np.zeros(
    shape=(in_df.index.shape[0], len(good_stns))).astype('float')
data_mtx[data_mtx == 0] = np.nan

df_all = pd.DataFrame(index=in_df.index, columns=good_stns, data=data_mtx)

# %%
df_stn0 = in_df.loc[:, good_stns].dropna(how='all')

#sorted_stn = stn0.sort_values()

#ata_index_sorted = np.squeeze(np.argsort(stn0.values, axis=0)[::-1])
#y0 = (np.arange(sorted_stn.size) / len(sorted_stn))

for stn_ in df_stn0.columns:

    # %% transform data to quantiles from cdf
    x01, y01 = build_edf_fr_vals(df_stn0.loc[:, stn_].values)
    df_stn_ = pd.DataFrame(data=df_stn0.loc[:, stn_].values,
                           index=df_stn0.loc[:, stn_].index, columns=[stn_])

    df_stn_.sort_values(by=stn_, inplace=True)

    df_stn_.loc[:, 'edf'] = y01
    df_stn_.sort_index(inplace=True)
    df_all.loc[:, stn_] = df_stn_.loc[:, 'edf']
    break
# %%
import matplotlib.pyplot as plt
# fig = plt.figure(figsize=(12, 8), dpi=300)
# plt.plot(x01, y01, c='r', alpha=.5)
# #plt.plot(x0, y0, c='b',  alpha=.5)
# plt.show()
#%%


# %%

#get_vals_fr_edf_vals(stn0, y0)
