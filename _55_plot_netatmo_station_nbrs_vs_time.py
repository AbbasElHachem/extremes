#-------------------------------------------------------------------------
# Name:        Resample
# Purpose:      resample precipatation data from natatmo weather stations
#
# Author:      Abbas El Hachem
#
# Created:
#
#-------------------------------------------------------------------------
#%%
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import time
import datetime as dt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.dates as md
# Set base folder
basefolder = r'X:\exchange\hespen\05_netatmo'
os.chdir(basefolder)

#%%
path_netatmo_ppt_data = r"ppt_all_netatmo_bw_hourly_everything.csv"

df = pd.read_csv(path_netatmo_ppt_data, sep=';', index_col=0,
                 parse_dates=True, infer_datetime_format=True)

df.dropna(how='all', inplace=True)

df = df.loc[:'2019-08-18 23:00:00', :]
#%%
#df.index = pd.to_datetime(df.index, format='%Y-%m-%d')

#%%
# resample from hour to month
#df_monthly = resampleDf(df, '1M')

# df2 = df.notnull().astype('int')
# sum_vals = df2.sum(axis=1)


# set all values >= 0 to 1
# df = df.apply(lambda x: [1 if y >= 0 else np.nan for y in x])
# df.fillna(1,)
# This is wrong because it sets np.nan equal to 1
#df = df.apply(lambda x: [y if y <= 0 else 1 for y in x])

#%%
# # sum up number of stations with data for each day
# df.loc[:, 'sum'] = pd.Series(df.sum(axis=1, skipna=True), index=df.index)
#
# # df_sum = df.loc[:, 'sum'][df.loc[:, 'sum'].values > 65]
# df_sum = df.loc[:, 'sum']
# df_sum = df_sum.iloc[:-1]


df_operation_times = pd.DataFrame(index=df.index)
dict_days = {ix: [0] for ix in df.index}
# find first valid idx per station
for i, stn in enumerate(df.columns):
    print(i, '/', len(df.columns))
    first_ix = df.loc[:, stn].notna().idxmax()
    dict_days[first_ix].append(1)


dict_days2 = {ix: [0] for ix in df.index}
for ix2 in dict_days2:
    dict_days2[ix2] = sum(dict_days[ix2])
df_days = pd.DataFrame.from_dict(dict_days2, orient='index')
vals = np.cumsum(df_days.values)


df3 = pd.DataFrame(index=df_days.index, data=vals)
# sum_vals2 = df3.sum(axis=1)


plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 14})

#%%
plt.ioff()
plt.figure(figsize=(12, 8), dpi=200)
ax = plt.subplot(111)

# ax.xaxis.set_major_formatter(myFmt)
# ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=10))
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()


xfmt = md.DateFormatter('%Y')

# ax.xaxis.set_major_locator(MultipleLocator(365))
# ax.yaxis.set_major_locator(MultipleLocator(2))

ax.xaxis.set_major_formatter(xfmt)
# ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.grid(color='k', linestyle='--', linewidth=0.1, alpha=0.5)

# ax.plot(['2015-12', '2015-12'],
#         [df3.min().values[0],
#          df3.max().values[0]], color='r', linestyle='--')

ax.set_ylabel('Netatmo Stations in BW with valid Observations')
# ax.set_xlabel('Time', ha="center")
ax.plot(df3.index, df3.values, color='b',
        alpha=0.75)

ax.axvline(x='2015-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
ax.axvline(x='2016-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
ax.axvline(x='2017-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
ax.axvline(x='2018-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
#ax.axvline(x='2019-12-31 00:00:00', color='k', linestyle='--', alpha=0.75)
ax.set_yticks(np.arange(0, 3200, 500))
# plt.xticks(rotation=45)
# plt.xticks([])
# Customize minor tick labels

ax.set_xticks(['2015-06-30 00:00:00',
               '2016-06-30 00:00:00',
               '2017-06-30 00:00:00',
               '2018-06-30 00:00:00',
               '2019-03-30 00:00:00'],
              minor=False)
# ax.set_xticklabels(['2015','2016','2017','2018','2019'], minor=True)

plt.xlim([df3.index[0] - pd.Timedelta(days=5),
          df3.index[-1] + pd.Timedelta(days=5)])
# plt.legend(loc='lower right')
# plt.tight_layout()
# can also save it as a PDF, JPEG, etc.
plt.savefig("netatmo_station_measurements2.png", papertype='a4',
            bbox_inches="tight")
plt.close('all')
