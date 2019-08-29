#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 08:08:15 2019

@author: abbas
"""

#%%
import pandas as pd
import numpy as np


#%%
df_netatmo_daily = pd.read_csv(r'/home/abbas/Desktop/data/data/all_netatmo_ppt_data_daily_.csv',
                               sep=';', index_col=0, parse_dates=True, infer_datetime_format=True,
                               engine='c').dropna(how='all')
#%%

find_maximum_dates = df_netatmo_daily.max(axis=1).sort_values()[::-1]

maximum_100_days = find_maximum_dates[:100]
#%%


#%%


#%%


#%%


#%%
