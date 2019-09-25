# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Look at pairs of stations (DWD-DWD) 

For different aggregations (15min, 30min, 60min, ..., 6hours, ..., 1 day)
Calculate P0 (probability that rainfall < 1mm)
Construct Cdf for extremes ( > 1mm) and compare
Calculate the correlation between the ranks on all scales
Find if the similarities or differences in the stations are due to 
a systematic or a random error.
If a systematic error is present (bias) find the correction factor 
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"
#==============================================================================
#
#==============================================================================
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from pandas.plotting import register_matplotlib_converters

from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from _00_additional_functions import (select_df_within_period,
                                      build_edf_fr_vals)

plt.ioff()

register_matplotlib_converters()

#==============================================================================
#
#==============================================================================
rc('font', size=13)
rc('font', family='serif')
rc('axes', labelsize=13)
rcParams['axes.labelpad'] = 13

majorLocator = MultipleLocator(2)
minorLocator = MultipleLocator(1)

path_to_ppt_dwd_data = (
    r"F:\download_DWD_data_recent\all_dwd_hourly_ppt_data_combined_1995_2019.fk")
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'


path_to_dwd_coords_only_in_bw = (
    r"F:\download_DWD_data_recent\station_coordinates_names_hourly_only_in_BW.csv")

out_save_dir_orig = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_DWD_compare_past_to_recent')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

x_col_name = 'Rechtswert'
y_col_name = 'Hochwert'

# threshold for CDF, consider only above thr, below is P0
ppt_thr_min = .5

# used for P0 calculation
ppt_thrs_list = [0.5]  # , 1, 2]


# till 1 day
aggregation_frequencies = ['60min', '120min', '480min',
                           '720min', '1440min']

old_period_start = '1995-01-01 00:00:00'
old_period_end = '2014-01-01 00:00:00'

new_period_start = '2014-01-01 00:00:00'
new_period_end = '2019-09-01 00:00:00'
#==============================================================================
#
#==============================================================================


def plot_end_tail_cdf_1_stn_old_recent(stn1_id, df1, df2, temp_freq,
                                       ppt_thr, out_dir):
    print('plotting CDF plots')

    values_x = df1.values.ravel()
    values_y = df2.values.ravel()

    start_df1, end_df1 = str(df1.index[0]), str(df1.index[-1])
    start_df2, end_df2 = str(df2.index[0]), str(df2.index[-1])

    # data_shape_df1 = np.round(df1.values.ravel().shape[0], 0)
    # data_shape_df2 = np.round(df2.values.ravel().shape[0], 0)

    xvals1, yvals1 = build_edf_fr_vals(values_x)
    xvals2, yvals2 = build_edf_fr_vals(values_y)

    fig = plt.figure(figsize=(20, 12), dpi=150)
    ax = fig.add_subplot(111)
    ax.scatter(xvals1, yvals1, c='b', marker='+', s=10,
               alpha=0.5,
               label='%s to %s' % (start_df1, end_df1))

    ax.scatter(xvals2, yvals2, c='r', marker='o', s=10,
               alpha=0.5,
               label='%s to %s' % (start_df2, end_df2))

    ax.set_xlim(min(xvals1.min(), xvals2.min()) - 0.1,
                max(xvals1.max(), xvals2.max()) + 1)
    #ax.set_yticks([0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1])
    ax.set_ylim(0.9, 1.005)

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    ax.set_xlabel('Observed %s Rainfall in mm station Id: %s'
                  % (temp_freq, stn1_id))
    ax.set_ylabel('Observed %s CDF station Id: %s'
                  % (temp_freq, stn1_id))
    ax.legend(loc='lower right')
    ax.set_title("Comapring CDF old and recent periods for DWD Stn: %s"
                 #"\n Available Data St1: %d; \n"
                 "\n Time Frequency: %s\n"
                 #"Available Data St2: %d; \n Time Frequency: %s\n"
                 "Rainfall Threshold: %.1f" % (stn1_id,
                                               # data_shape_df1,
                                               # data_shape_df2,
                                               temp_freq,
                                               ppt_thr))

    ax.grid(color='k', linestyle='--', linewidth=0.5, alpha=0.25)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,
                             'compare_old_vs_recent_%s_cdf_stn_%s_.png'
                             % (temp_freq, stn1_id)))
    plt.clf()
    plt.close('all')
    print('Done saving figure CDF')

    return


#==============================================================================
#
#==============================================================================
print('Reading Coordinate Dataframes')

in_coords_df = pd.read_csv(path_to_dwd_coords_only_in_bw,
                           sep=',', index_col=0, encoding='latin')

stndwd_ix = ['0' * (5 - len(str(stn_id))) + str(stn_id)
             if len(str(stn_id)) < 5 else str(stn_id)
             for stn_id in in_coords_df.index]

in_coords_df.index = stndwd_ix

print('Reading Data Dataframes')

df_dwd = pd.read_feather(path_to_ppt_dwd_data, use_threads=True)
df_dwd.set_index('Time', inplace=True)

df_dwd.index = pd.to_datetime(df_dwd.index, format='%Y-%m-%d %H:%M:%S')
df_dwd.dropna(axis=0, inplace=True, how='all')

print('Splitting Data Dataframes')
df_dwd_old = select_df_within_period(df=df_dwd,
                                     start=old_period_start,
                                     end=old_period_end)

df_dwd_new = select_df_within_period(df=df_dwd,
                                     start=new_period_start,
                                     end=new_period_end)
print('Start plotting')

for dwd_stn_id in df_dwd_old.columns:
    if dwd_stn_id in stndwd_ix:
        print('plotting for', dwd_stn_id)
        df_old = df_dwd_old.loc[:, dwd_stn_id].dropna(how='all')
        df_new = df_dwd_new.loc[:, dwd_stn_id].dropna(how='all')
        if len(df_old.values) > 0 and len(df_new.values) > 0:
            plot_end_tail_cdf_1_stn_old_recent(stn1_id=dwd_stn_id,
                                               df1=df_old,
                                               df2=df_new,
                                               temp_freq='60min',
                                               ppt_thr=0,
                                               out_dir=out_save_dir_orig)

#            break
