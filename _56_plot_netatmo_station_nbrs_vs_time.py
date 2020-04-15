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


plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'axes.labelsize': 26})


from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

# for BW
# Set base folder
# basefolder = r'X:\exchange\hespen\05_netatmo'
# os.chdir(basefolder)
#
# #%%
# path_netatmo_ppt_data = r"ppt_all_netatmo_bw_hourly_everything.csv"
#
# path_netatmo_ppt_data_daily = (r"X:\hiwi\ElHachem\Prof_Bardossy"
#                                r"\Extremes\NetAtmo_BW\ppt_all_netatmo_1440min_.csv")

# for RLP
basefolder = r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz'
os.chdir(basefolder)

#%%
path_netatmo_ppt_data = r"ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.csv"


df = pd.read_csv(path_netatmo_ppt_data, sep=';', index_col=0, engine='c',
                 parse_dates=True, infer_datetime_format=True)

df.dropna(how='all', inplace=True)

df = df.loc['2014-06-01 00:00:00':'2019-10-31 23:00:00', :]
# '2019-08-31 23:00:00'

#%%
#df.index = pd.to_datetime(df.index, format='%Y-%m-%d')

#==============================================================================
# observation wit time
#==============================================================================
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

#==============================================================================
# observation per station
#==============================================================================
# sum_vals2 = df3.sum(axis=1)
df_netatmo_ppt = pd.read_csv(path_netatmo_ppt_data,
                             # path_netatmo_ppt_data_daily, # for BW
                             sep=';', index_col=0,
                             parse_dates=True, infer_datetime_format=True)
df_netatmo_ppt.dropna(how='all', inplace=True)


nbr_days_per_stn = []
for stn in df_netatmo_ppt.columns:
    df_stn = df_netatmo_ppt.loc[:, stn].dropna(how='all')
    nbr_days_per_stn.append(df_stn.size / 24)  # for BW not / 24

nbr_days_per_stn = np.array(nbr_days_per_stn)

number_of_data_per_stn = nbr_days_per_stn

df_count = pd.Series(data=number_of_data_per_stn,
                     index=range(1, number_of_data_per_stn.shape[0] + 1))

# df_count.plot(legend=False)

#==============================================================================
# PLOT
#==============================================================================
plt.ioff()
fig, axs = plt.subplots(1, 2, sharex=False, sharey=False,
                        figsize=(32, 16), dpi=300)
# (31, 12)

bins = [0, 360, 720, 1095, 1460, 1825]

hist_, bin_edges = np.histogram(df_count,
                                bins=bins, density=False)
hist_vals = [v for v in hist_]
axs[1].bar(bins[1:], hist_vals, facecolor='b', alpha=0.95, width=50,
           align='center',
           ec="k", lw=1.)
#
axs[1].set_yticks([0, 300, 600, 900])

#     np.arange(0, 1100, 200))
# axs[1].set_xlim([0, 1830])

#axs[1].set_ylabel("Number of Netatmo stations")
axs[1].set_xlabel("Years with valid observations", labelpad=16)

# axs[1].axvline(x=180, color='k', linestyle='--', alpha=0.05)
axs[1].axvline(x=540, color='k', linestyle='--', alpha=0.1)
axs[1].axvline(x=907.5, color='k', linestyle='--', alpha=0.1)
axs[1].axvline(x=1277.5, color='k', linestyle='--', alpha=0.1)
axs[1].axvline(x=1642.5, color='k', linestyle='--', alpha=0.1)
axs[1].set_xticks(bins[1:])

#[180, 540, 907.5, 1277.5, 1642.5]
# axs[0].set_xticklabels(['< 1 [y]', '1-2 [y]', '2-3 [y]', '3-4 [y]', '>4 [y]'])
axs[1].set_xticklabels(['< 1', '1-2', '2-3', '3-4', '>4'])


axs[1].grid(color='k', linestyle='--', linewidth=0.1, alpha=0.5)

# ax.plot(['2015-12', '2015-12'],
#         [df3.min().values[0],
#          df3.max().values[0]], color='r', linestyle='--')

axs[0].set_ylabel(
    'Netatmo Stations in RLP with valid Observations', labelpad=16)  # BW
# ax.set_xlabel('Time', ha="center")
axs[0].plot(df3.index, df3.values, color='b',
            alpha=0.95)
xfmt = md.DateFormatter('%Y')

axs[0].xaxis.set_major_formatter(xfmt)
# ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

axs[0].grid(color='k', linestyle='--', linewidth=0.1, alpha=0.5)
axs[0].axvline(x='2014-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
axs[0].axvline(x='2015-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
axs[0].axvline(x='2016-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
axs[0].axvline(x='2017-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
axs[0].axvline(x='2018-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)
axs[0].axvline(x='2019-12-31 00:00:00', color='k', linestyle='--', alpha=0.25)


#ax.axvline(x='2019-12-31 00:00:00', color='k', linestyle='--', alpha=0.75)
axs[0].set_yticks([0, 500, 1000, 1500, 2000, 2500])  # 3000
# np.arange(0, 3200, 250))

axs[0].set_xticks(['2014-06-30 00:00:00',
                   '2015-06-30 00:00:00',
                   '2016-06-30 00:00:00',
                   '2017-06-30 00:00:00',
                   '2018-06-30 00:00:00',
                   '2019-06-30 00:00:00'],  # 03 bor BW
                  minor=False)
axs[0].set_xlim([df3.index[0] - pd.Timedelta(days=5),
                 df3.index[-1] + pd.Timedelta(days=5)])
# ax.set_xticklabels(['2015','2016','2017','2018','2019'], minor=True)

axs[0].legend(title='a)', loc='upper left',  # a) left
              frameon=False, fontsize=26)._legend_box.align = 'left'

axs[1].legend(title='b)', loc='upper left',
              frameon=False, fontsize=26)._legend_box.align = 'left'


# plt.legend(loc='lower right')
# for BW
# X:\hiwi\ElHachem\Prof_Bardossy\Extremes\netatmo_station_measurements_2.png"
plt.tight_layout()
# can also save it as a PDF, JPEG, etc.
plt.savefig(
    r"RLP_data.png",
    papertype='a4',
    bbox_inches="tight")
plt.close('all')
