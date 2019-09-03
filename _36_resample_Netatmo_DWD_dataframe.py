
# coding: utf-8

# In[13]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

#path_to_df = r"F:\download_DWD_data_recent\all_dwd_hourly_ppt_data_combined_1995_2019.csv"
path_to_df = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new.csv"
#path_to_df = r"F:\download_DWD_data_recent\all_dwd_10min_ppt_data_combined_2014_2019.csv"


# In[2]:


def agg_data(df, agg, closed, label, shift=False, leave_nan=True,
             max_nan=0):
    """Aggregate precipitation data

    Parameters:
    -----------
    Df: Pandas DataFrame Object
        Data set
    agg: string
        Aggregation 'M':Monthly 'D': Daily, 'H': Hourly, 'Min': Minutely
    closed: string
        'left' or 'right' defines the aggregation interval
    label: string
        'left' or 'right' defines the related timestamp
    shift: boolean, optional
        Shift the values by 6 hours according to the dwd daily station.
        Only valid for aggregation into daily aggregations
        True, data is aggregated from 06:00 - 06:00
        False, data is aggregated from 00:00 - 00:00
        Default is False

    leave_nan: boolean, optional
        True, if the nan values should remain in the aggregated data set.
        False, if the nan values should be treated as zero values for the
        aggregation. Default is True

    Remark:
    -------
        If the timestamp is at the end of the timeperiod:

        Input: daily        Output: daily+
            >> closed='right', label='right'

        Input: subdaily     Output: subdaily
            >> closed='right', label='right'

        Input: subdaily     Output: daily
            >> closed='right', label='left'

        Input: subdaily     Output: monthly
            >> closed='right', label='right'


        ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
        ! ! Always check, if aggregation is correct ! !
        ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !



    """

    if shift == True:
        df_copy = df.copy()
        if agg != 'D':
            raise Exception('Shift can only be applied to daily aggregations')
        df = df.shift(-6, 'H')

    # To respect the nan values
    if leave_nan == True:
        # for max_nan == 0, the code runs faster if implemented as follows
        if max_nan == 0:
            # Fill the nan values with values very great negative values and later
            # get the out again, if the sum is still negative
            df = df.fillna(-100000000000.)
            df_agg = df.resample(agg,
                                 closed=closed,
                                 label=label).sum()
            # Replace negative values with nan values
            df_agg.values[df_agg.values[:] < 0.] = np.nan
        else:
            df_agg = df.resample(rule=agg,
                                 closed=closed,
                                 label=label).sum()
            # find data with nan in original aggregation
            g_agg = df.groupby(pd.Grouper(freq=agg,
                                          closed=closed,
                                          label=label))
            n_nan_agg = g_agg.aggregate(lambda x: pd.isnull(x).sum())

            # set aggregated data to nan if more than max_nan values occur in the
            # data to be aggregated
            filter_nan = (n_nan_agg > max_nan)
            df_agg[filter_nan] = np.nan

    elif leave_nan == False:
        df_agg = df.resample(agg,
                             closed=closed,
                             label=label).sum()
    if shift == True:
        df = df_copy
    return df_agg


# In[3]:


df_hourly = pd.read_csv(path_to_df, sep=';', index_col=0,
                        infer_datetime_format=True, engine='c', parse_dates=True).dropna(how='all')


# In[4]:


df_daily = agg_data(df_hourly, '60min', closed='right',
                    label='right', max_nan=0)


# In[9]:


df_daily = df_daily[0 <= df_daily]
df_daily = df_daily[df_daily < 300]
df_daily.dropna(how='all', inplace=True)


# In[20]:


# all_netatmo_ppt_data_daily_
# all_dwd_hourly_ppt_data_combined_2014_2019_
df_daily.to_csv(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\all_netatmo_ppt_data_daily_.csv",
                sep=';', float_format='%.3f')
