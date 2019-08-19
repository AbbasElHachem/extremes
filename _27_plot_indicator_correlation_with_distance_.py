# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Plot Spearman Indicator Correlation with Distance
Purpose: Filter data based on spearman indicator correlation

Created on: 2019-08-10

After calculating for every Netatmo station the Spearman indicator rank 
correlation with either neigboring DWD Netatmo station (neighbors 1 to 5)
plot the result, x axis for distance and y axis for correlation
find the distribution of correaltion with increasing distance and
choice of neighbor.

This is used as a filter to remove the 'bad' Netamo staions 


Parameters
----------

Input Files
    DF file for each neighbor, index is station name, columns contain 
    distance to neighbor and rank spearman correlation
    
    Either Netatmo-DWD or Netatmo-Netatmo or DWD-DWD 
    for different percentage threshold and different time aggregations
Returns
-------

    Plot correlation with distance for different neighbors
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')

data_dir_Netamto_dfs = main_dir / r'plots_NetAtmo_ppt_DWD_ppt_correlation_'

assert data_dir_Netamto_dfs.exists(), 'Wrong Netatmo path'


data_dir_DWD_dfs = main_dir / r'plots_DWD_ppt_DWD_ppt_correlation_'
assert data_dir_DWD_dfs.exists(), 'Wrong dwd path'

data_dir_Netamto_netatmo_dfs = main_dir / \
    r'plots_NetAtmo_ppt_Netatmo_ppt_correlation_'

assert data_dir_Netamto_netatmo_dfs.exists(), 'Wrong Netatmo Netatmo path'

netatmo_path_acc = r'year_allyears_df_comparing_correlations_max_sep_dist_500000_'
dwd_path_Acc = r'year_allyears_df_dwd_correlations'

# def percentage threshold, time frequency and data source
percent = '85'
time_freq = '1440min'

data_source0 = 'DWD'  # 'Netatmo'  #   # reference station 'Netatmo'
data_source = 'dwd'  # 'netatmo'  #   # compare to station 'netatmo'

# =============================================================================


def gen_path_df_file(main_path, start_path_acc, time_freq,
                     data_source, percent, neighbor):
    ''' use this funtion to get the path to the dfs for different neighbors'''

    return main_path / (
        r'%sfreq_%s_%s_netatmo_upper_%s_percent_data_considered_neighbor_%d_.csv'
        % (start_path_acc, time_freq, data_source, percent, neighbor))


# =============================================================================
if data_source0 == 'Netatmo' and data_source == 'dwd':
    # for Netatmo stations neighbors start from 0 !
    df0 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 0)
    df1 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 1)
    df2 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 2)
    df3 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 3)
    df4 = gen_path_df_file(data_dir_Netamto_dfs, netatmo_path_acc, time_freq,
                           data_source, percent, 4)
    save_dir = data_dir_Netamto_dfs

if data_source0 == 'DWD' and data_source == 'dwd':
    # for DWD stations neighbors start from 1 not 0  (1 is first)!
    df0 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 1)
    df1 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 2)
    df2 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 3)
    df3 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 4)
    df4 = gen_path_df_file(data_dir_DWD_dfs, dwd_path_Acc, time_freq,
                           data_source, percent, 5)
    save_dir = data_dir_DWD_dfs

if data_source0 == 'Netatmo' and data_source == 'netatmo':
    # for Netatmo stations neighbors start from 0 !
    df0 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
                           netatmo_path_acc, time_freq,
                           data_source, percent, 0)
    df1 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
                           netatmo_path_acc, time_freq,
                           data_source, percent, 1)
    df2 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
                           netatmo_path_acc, time_freq,
                           data_source, percent, 2)
    df3 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
                           netatmo_path_acc, time_freq,
                           data_source, percent, 3)
    df4 = gen_path_df_file(data_dir_Netamto_netatmo_dfs,
                           netatmo_path_acc, time_freq,
                           data_source, percent, 4)
    save_dir = data_dir_Netamto_netatmo_dfs
# =============================================================================

in_df0 = pd.read_csv(df0, index_col=0, sep=';').dropna(how='all')
in_df1 = pd.read_csv(df1, index_col=0, sep=';').dropna(how='all')
in_df2 = pd.read_csv(df2, index_col=0, sep=';').dropna(how='all')
in_df3 = pd.read_csv(df3, index_col=0, sep=';').dropna(how='all')
in_df4 = pd.read_csv(df4, index_col=0, sep=';').dropna(how='all')

# =============================================================================

x0 = in_df0.loc[:, 'Distance to neighbor'].values.ravel()
x1 = in_df1.loc[:, 'Distance to neighbor'].values.ravel()
x2 = in_df2.loc[:, 'Distance to neighbor'].values.ravel()
x3 = in_df3.loc[:, 'Distance to neighbor'].values.ravel()
x4 = in_df4.loc[:, 'Distance to neighbor'].values.ravel()

y0 = in_df0.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y1 = in_df1.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y2 = in_df2.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y3 = in_df3.loc[:, 'Bool_Spearman_Correlation'].values.ravel()
y4 = in_df4.loc[:, 'Bool_Spearman_Correlation'].values.ravel()

# =============================================================================

plt.ioff()
plt.figure(figsize=(12, 8), dpi=300)
plt.scatter(x0, y0, c='r', alpha=0.5, marker='x',
            label='First Neighbor Stn nbr %d' % y0.shape[0], s=28)
plt.scatter(x1, y1, c='b', alpha=0.5, marker='.',
            label='Second Neighbor Stn nbr %d' % y1.shape[0], s=28)
plt.scatter(x2, y2, c='g', alpha=0.5, marker='d',
            label='Third Neighbor Stn nbr %d' % y2.shape[0], s=28)
plt.scatter(x3, y3, c='orange', alpha=0.5, marker='1',
            label='Fourth Neighbor Stn nbr %d' % y3.shape[0], s=28)
plt.scatter(x4, y4, c='m', alpha=0.5, marker='2',
            label='Fifth Neighbor Stn nbr %d' % y4.shape[0], s=28)

plt.xlim([0, max(x3.max(), x4.max()) + 1000])
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
plt.tight_layout()
plt.title('%s %s stations %s: Indicator correlation'
          ' with distance for upper %s percent of data values'
          % (data_source0, data_source, time_freq, percent))
plt.savefig(save_dir /
            (r'_%s_%s_%s_percent_indic_corr_freq_%s_.png'
             % (data_source0, data_source, percent, time_freq)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
plt.close()
# plt.show()
