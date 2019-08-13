#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:58:27 2019

@author: abbas
"""
#%%
import pandas as pd
import numpy as np
# create indicator Kriging

df_file = r'/home/abbas/Desktop/home-office/ppt_all_netatmo_hourly_stns_combined_.fk'
#%%
in_df = pd.read_feather(df_file, use_threads=True)
in_df.set_index('Time', inplace=True)
in_df.index = pd.to_datetime(in_df.index, format='%Y-%m-%d %H:%M:%S')
#%%

all_stns = in_df.columns.to_list()
# %%
stn0 = all_stns[100]
df_stn0 = in_df.loc[:, stn0].dropna()

#sorted_stn = stn0.sort_values()

#ata_index_sorted = np.squeeze(np.argsort(stn0.values, axis=0)[::-1])
#y0 = (np.arange(sorted_stn.size) / len(sorted_stn))

# %% transform data to quantiles from cdf



def build_edf_fr_vals(ppt_data):
    # Construct EDF
    ''' construct empirical distribution function given data values '''
    data_sorted = np.sort(ppt_data, axis=0)[::-1]
    data_index_sorted = np.squeeze(np.argsort(ppt_data, axis=0)[::-1])
    x0 = np.squeeze(data_sorted)[::-1]
    y0 = (np.arange(data_sorted.size) / len(data_sorted))
    return x0, y0

x0, y0 = build_edf_fr_vals(df_stn0.values)

# %%
def get_edf2(data):
    # create a sorted series of unique data
    # cdfx = np.sort(np.unique(data))
    # x-data for the ECDF: evenly spaced sequence of the uniques
    data_sorted = np.sort(data, axis=0)
    data_index_sorted = np.squeeze(np.argsort(data, axis=0))
    # size of the x_values
    size_data = data.size
    # y-data for the ECDF:
    y_values = []

    for i in data_sorted:
        # all the values in raw data less than the ith value in x_values
        temp = data[data <= i]
        # fraction of that value with respect to the size of the x_values
        value = np.round(temp.size / size_data, 3)
        # pushing the value in the y_values
        y_values.append(value)

    # return both x and y values    
    return data_sorted, y_values, data_index_sorted

x01, y01, ix = get_edf2(df_stn0.values)

# %%
df_stn0_ = pd.DataFrame(index=df_stn0.index, data=df_stn0.values, columns=[stn0])
df_stn0_.sort_values(by=stn0, inplace=True)
df_stn0_.loc[:,'edf'] = y01
df_stn0_.sort_index(inplace=True)

# %%
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(12,8), dpi=300)
plt.plot(x01, y01, c='r', alpha=.5)
#plt.plot(x0, y0, c='b',  alpha=.5)
plt.show()
#%%


# %%

#get_vals_fr_edf_vals(stn0, y0)