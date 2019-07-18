# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Check Tendency of a Station
Purpose: See if a Netatmo station has same behavior in p1, p5, mean_vals,
     correlations as compared to it's neighbouring DWD station

Created on: 2019-07-18

Parameters
----------
file_loc : str
    The file location of the spreadsheet
print_cols : bool, optional
    A flag used to print the columns to the console (default is False)

Returns
-------
list
    a list of strings representing the header columns
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================

from pathlib import Path

import os
import timeit
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from _00_additional_functions import list_all_full_path


main_dir = Path(os.getcwd())
os.chdir(main_dir)

path_to_netatmo_df = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                      '\plots_NetAtmo_ppt_NetAtmo_temperature')


ppt_thr = 1
dfs = list_all_full_path('.csv', path_to_netatmo_df)

dfs_tendency_p1_p5 = [df_file for df_file in dfs
                      if ('p%d_and_mean' % ppt_thr) in df_file]

dfs_tendency_correlations = [df_file for df_file in dfs
                             if ('correlations_p%d' % ppt_thr) in df_file]

stations_all = []
for df_file in dfs_tendency_p1_p5:
    df_c = pd.read_csv(df_file, sep=';', index_col=0)
    for index, row in df_c.iterrows():
        stations_all.append(index)
        #print(row['p1'], row['mean_ppt'])
        #print(row['Orig_Pearson_Correlation'], row['Orig_Spearman_Correlation'] , Bool_Pearson_Correlation)


stations_all_uniq = np.unique(stations_all).astype(str)
dfs_years = pd.DataFrame(index=stations_all_uniq,
                         columns=['all_years', '2015', '2016', '2017', '2018', '2019'])

for year in dfs_years.columns:
    for df_f in dfs_tendency_p1_p5:
        if year in df_f:
            print(year, df_f)
            df_c = pd.read_csv(df_f, sep=';', index_col=0)
            for index, row in df_c.iterrows():
                if index not in dfs_years.index:
                    print(df_c.loc['70_ee_50_00_42_e2', :])
                dfs_years.loc[index, year] = row['Orig_Pearson_Correlation']

dfs_year_num = dfs_years.replace('_', -1).replace('s', 0).replace('+', 1)
dfs_year_num.to_csv(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\dfs_comparing\df_tendency_Orig_Pearson_Correlation.csv', sep=';')


plt.ioff()
fig = plt.figure(figsize=(16, 8), dpi=200)

ax = fig.add_subplot(111)

markes = ['+', 'o', ',', '<', '1']
colors = ['r', 'b', 'g', 'k', 'm']

# for i, year in enumerate(dfs_year_num.columns):
#    idx = dfs_year_num.loc[:, year].index
#    vals = dfs_year_num.loc[:, year].values
#    ax.scatter(idx, vals, c=colors[i], alpha=0.15,
#               marker=markes[i], s=5, label=year)
# break

for i, station in enumerate(dfs_year_num.index):
    idx = dfs_year_num.columns
    vals = dfs_year_num.loc[station, :].values

    ax.plot(idx, vals, c='k', alpha=0.5,
            marker=',', linewidth=1, label=station)

plt.ylim([0, 1.1])

# plt.axis('equal')
plt.grid(alpha=.005)
#plt.legend(loc=0, fontsize=12)
plt.xlabel('Netatmo Station', fontsize=12)
plt.ylabel('Changes in Pearson Correlation', fontsize=12)
plt.xticks(fontsize=12, rotation=15)
plt.yticks([0, 1], fontsize=12)
# pt.tight_layout()
plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\corr_along_years.png",
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)

if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
