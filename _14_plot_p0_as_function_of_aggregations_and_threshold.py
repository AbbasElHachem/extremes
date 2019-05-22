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
import timeit
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy.spatial as spatial


from matplotlib import rc
from matplotlib import rcParams

from _00_additional_functions import (resampleDf,
                                      get_cdf_part_abv_thr)

from b_get_data import HDF5


plt.ioff()
#==============================================================================
#
#==============================================================================
rc('font', size=13)
rc('font', family='serif')
rc('axes', labelsize=13)
rcParams['axes.labelpad'] = 13


path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')

coords_df_file = (r'X:\exchange\ElHachem\niederschlag_deutschland'
                  r'\1993_2016_5min_merge_nan.csv')
assert os.path.exists(coords_df_file), 'wrong DWD coords file'

out_save_dir = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_DWD_DWD')


if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)

x_col_name = 'Rechtswert'
y_col_name = 'Hochwert'

# threshold for CDF, consider only above thr, below is P0
ppt_thrs_list = [0.25, 0.5, 1, 2]

# till 1 day
aggregation_frequencies = ['5min', '10min', '15min', '30min', '60min', '90min',
                           '120min', '180min', '240min',  '360min',
                           '480min', '720min', '1440min']


#==============================================================================
#
#==============================================================================


def get_dwd_stns_coords(coords_df_file, x_col_name, y_col_name):
    '''function used to return to coordinates dataframe of the DWD stations'''
    in_coords_df = pd.read_csv(coords_df_file, sep=';',
                               index_col=3, engine='c')
    stn_ids = in_coords_df.index
    x_vals = in_coords_df[x_col_name].values.ravel()
    y_vals = in_coords_df[y_col_name].values.ravel()
    return in_coords_df, x_vals, y_vals, stn_ids
#==============================================================================
#
#==============================================================================


def get_nearest_dwd_station(first_stn_id, coords_df_file,
                            x_col_name, y_col_name):
    ''' Find for one station, the closest neibhouring station'''
    # read df coordinates and get station ids, and x, y values
    in_coords_df, x_vals, y_vals, stn_ids = get_dwd_stns_coords(
        coords_df_file, x_col_name, y_col_name)
    # make a tupples from the coordinates
    coords_tuples = np.array([(x, y) for x, y in zip(x_vals, y_vals)])
    # create a tree from coordinates
    points_tree = spatial.cKDTree(coords_tuples)

    xstn = in_coords_df.loc[int(first_stn_id), x_col_name]
    ystn = in_coords_df.loc[int(first_stn_id), y_col_name]

    distances, indices = points_tree.query([xstn, ystn], k=2)
    coords_nearest_nbr = coords_tuples[indices[1]]
    stn_near = str(stn_ids[indices[1]])
    distance_near = distances[1]

    return coords_nearest_nbr, stn_near, distance_near
#==============================================================================
#
#==============================================================================


def get_p0_for_2_stns_temp_freq(df1, df2, ppt_thr):
    print('Getting Probabilities P01, P02')
    values_x = df1.values.ravel()
    values_y = df2.values.ravel()

    _, yvals1 = get_cdf_part_abv_thr(values_x, ppt_thr)
    _, yvals2 = get_cdf_part_abv_thr(values_y, ppt_thr)

    p01 = yvals1[0]
    p02 = yvals2[0]

    return p01, p02
#==============================================================================
#
#==============================================================================


def compare_cdf_two_dwd_stns(stns_ids):
    for iid in stns_ids[10:]:
        print('First Stn Id is', iid)
        try:
            idf1 = HDF52.get_pandas_dataframe(ids=[iid])
            idf1.dropna(axis=0, inplace=True)

            _, stn_near, distance_near = get_nearest_dwd_station(
                iid, coords_df_file, x_col_name, y_col_name)
            assert iid != stn_near, 'wrong neighbour selected'
            idf2 = HDF52.get_pandas_dataframe(ids=stn_near)
            print('Second Stn Id is', stn_near)
            df_p01_stn1 = pd.DataFrame(columns=aggregation_frequencies,
                                       index=ppt_thrs_list)
            df_p01_stn2 = pd.DataFrame(columns=aggregation_frequencies,
                                       index=ppt_thrs_list)
            for ppt_thr in ppt_thrs_list:
                print('Ppt Threshold is', ppt_thr)
                for tem_freq in aggregation_frequencies:
                    print('Aggregation is: ', tem_freq)
                    df_resample1 = resampleDf(data_frame=idf1,
                                              temp_freq=tem_freq)
                    df_resample2 = resampleDf(data_frame=idf2,
                                              temp_freq=tem_freq)

                    idx_common = df_resample1.index.intersection(
                        df_resample2.index)

                    df_common1 = df_resample1.loc[idx_common, :]
                    df_common2 = df_resample2.loc[idx_common, :]
                    print('done interscting Dataframes')
                    if (df_common1.values.shape[0] > 0 and
                            df_common2.values.shape[0] > 0):
                        try:
                            p01, p02 = get_p0_for_2_stns_temp_freq(df_common1,
                                                                   df_common2,
                                                                   ppt_thr)
                            df_p01_stn1.loc[ppt_thr, tem_freq] = p01
                            df_p01_stn2.loc[ppt_thr, tem_freq] = p02

                        except Exception as msg:
                            print('error while plotting', msg, tem_freq)
                            continue
                            # break
                    else:
                        print('empty df')
                        continue
                print('plotting for Ppt thr', ppt_thr)

                plt.plot(aggregation_frequencies,
                         df_p01_stn1.loc[ppt_thr, :], c='r',
                         marker='o', linestyle='--',
                         label=ppt_thr, alpha=0.5)
                plt.plot(aggregation_frequencies,
                         df_p01_stn2.loc[ppt_thr, :], c='b',
                         marker='+', linestyle='-.', alpha=0.5,
                         label=ppt_thr)
            plt.legend(loc=0)
            plt.show()
        except Exception as msg:
            print(msg)

        break


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    HDF52 = HDF5(infile=path_to_ppt_hdf_data)
    ids = HDF52.get_all_ids()
    compare_cdf_two_dwd_stns(ids)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
