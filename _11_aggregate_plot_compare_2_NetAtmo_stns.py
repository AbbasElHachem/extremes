# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Look at pairs of stations (NetAtmo-NetAtmo) 

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

import scipy.spatial as spatial

from _00_additional_functions import (resample_intersect_2_dfs)

from _10_aggregate_plot_compare_2_DWD_stns import (plt_bar_plot_2_stns,
                                                   plt_scatter_plot_2_stns,
                                                   plot_end_tail_cdf_2_stns,
                                                   plot_normalized_ranked_stns)


path_to_ppt_netatmo_data = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                            r'\NetAtmo_BW'
                            r'\ppt_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_ppt_netatmo_data), 'wrong NETATMO Ppt file'

coords_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                  r'\NetAtmo_BW\rain_bw_1hour\netatmo_bw_1hour_coords.csv')

assert os.path.exists(coords_df_file), 'wrong NETATMO coords file'

out_save_dir = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_NetAtmo_Netatmo')


if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)

x_col_name = ' lon'
y_col_name = ' lat'


# threshold for CDF, consider only above thr, below is P0
ppt_thr = .5
max_ppt_thr = 100.

ppt_thr_2 = 0.5  # for the scatter plot

# till 1 day '5min', '10min', '15min', '30min',
aggregation_frequencies = ['60min',
                           '90min', '120min', '180min', '240min',
                           '360min', '480min', '720min', '1440min']
#==============================================================================
#
#==============================================================================


def get_netatmo_stns_coords(coords_df_file, x_col_name, y_col_name):
    '''function used to return to coordinates dataframe of the DWD stations'''
    in_coords_df = pd.read_csv(coords_df_file, sep=';',
                               index_col=0, engine='c')
    stn_ids = list(set([stn_id.replace(':', '_')
                        for stn_id in in_coords_df.index]))

    lon_vals = in_coords_df[x_col_name].values.ravel()
    lat_vals = in_coords_df[y_col_name].values.ravel()

    return in_coords_df, lon_vals, lat_vals, stn_ids
#==============================================================================
#
#==============================================================================


def get_nearest_netatmo_station(first_stn_id, coords_df_file,
                                x_col_name, y_col_name):
    ''' Find for one station, the closest neibhouring station'''
    # read df coordinates and get station ids, and x, y values
    in_coords_df, x_vals, y_vals, stn_ids = get_netatmo_stns_coords(
        coords_df_file, x_col_name, y_col_name)

    first_stn_id = first_stn_id.replace('_', ':')
    lon_stn = in_coords_df.loc[first_stn_id, x_col_name]
    lat_stn = in_coords_df.loc[first_stn_id, y_col_name]

    # make a tupples from the coordinates
    coords_tuples = np.array([(x, y) for x, y in zip(x_vals, y_vals)])

    # create a tree from coordinates
    points_tree = spatial.cKDTree(coords_tuples)

    distances, indices = points_tree.query([lon_stn, lat_stn], k=2)

    coords_nearest_nbr = coords_tuples[indices[1]]
    stn_near = str(stn_ids[indices[1]])
    distance_near = distances[1]

    return coords_nearest_nbr, stn_near, distance_near

#==============================================================================
#
#==============================================================================


def compare_cdf_two_stns(netatmo_ppt_df_file):
    in_netatmo_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                     index_col=0, sep=';',
                                     parse_dates=True,
                                     infer_datetime_format=True,
                                     engine='c')
    stns_ids = in_netatmo_stns_df.columns

    for stn_id in stns_ids[5:]:
        print('First Stn Id is', stn_id)
        try:
            idf1 = in_netatmo_stns_df.loc[:, stn_id]
            idf1.dropna(axis=0, inplace=True)
            idf1 = idf1[idf1 < max_ppt_thr]
            _, stn_near, distance_near = get_nearest_netatmo_station(
                stn_id, coords_df_file, x_col_name, y_col_name)
            assert stn_id != stn_near, 'wrong neighbour selected'
            idf2 = in_netatmo_stns_df.loc[:, stn_near]
            idf2 = idf2[idf2 < max_ppt_thr]
            print('Second Stn Id is', stn_near)
            for tem_freq in aggregation_frequencies:
                print('Aggregation is: ', tem_freq)
                df_common1, df_common2 = resample_intersect_2_dfs(idf1,
                                                                  idf2,
                                                                  tem_freq)
                if (df_common1.values.shape[0] > 10 and
                        df_common2.values.shape[0] > 10):
                    try:
                        plt_bar_plot_2_stns(stn_id,
                                            stn_near,
                                            distance_near,
                                            df_common1,
                                            df_common2,
                                            tem_freq,
                                            out_save_dir)
                        plt_scatter_plot_2_stns(stn_id,
                                                stn_near,
                                                distance_near,
                                                df_common1,
                                                df_common2,
                                                ppt_thr_2,
                                                tem_freq,
                                                out_save_dir)
                        plot_end_tail_cdf_2_stns(stn_id,
                                                 stn_near,
                                                 distance_near,
                                                 df_common1,
                                                 df_common2,
                                                 tem_freq,
                                                 ppt_thr,
                                                 out_save_dir)
                        plot_normalized_ranked_stns(stn_id,
                                                    stn_near,
                                                    distance_near,
                                                    df_common1,
                                                    df_common2,
                                                    tem_freq,
                                                    out_save_dir)
                    except Exception as msg:
                        print('error while plotting', msg, tem_freq)
                        continue
                    # break
                else:
                    print('empty df')
                    continue
        except Exception as msg:
            print(msg)

        break


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    compare_cdf_two_stns(path_to_ppt_netatmo_data)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
