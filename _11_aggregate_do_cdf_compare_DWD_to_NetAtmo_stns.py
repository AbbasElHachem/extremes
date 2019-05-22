# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Look at pairs of stations (DWD-NetAtmo) 

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


from _00_additional_functions import (
    convert_coords_fr_wgs84_to_utm32_, resampleDf)


from _09_aggregate_do_cdf_compare_2_DWD_stns import (plt_bar_plot_2_stns,
                                                     plt_scatter_plot_2_stns,
                                                     plot_end_tail_cdf_2_stns,
                                                     plot_ranked_stns)

from b_get_data import HDF5
#==============================================================================
#
#==============================================================================

path_to_ppt_netatmo_data = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                            r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_ppt_netatmo_data), 'wrong NETATMO Ppt file'

path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')
assert os.path.exists(path_to_ppt_hdf_data), 'wrong NETATMO Ppt file'

coords_dwd_df_file = (r'X:\exchange\ElHachem\niederschlag_deutschland'
                      r'\1993_2016_5min_merge_nan.csv')
assert os.path.exists(coords_dwd_df_file), 'wrong DWD coords file'

coords_netatmo_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                          r'\NetAtmo_BW\rain_bw_1hour\netatmo_bw_1hour_coords.csv')

assert os.path.exists(coords_netatmo_df_file), 'wrong NETATMO coords file'

out_save_dir = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_DWD_NetAtmo')


if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)

lon_col_name = ' lon'
lat_col_name = ' lat'

x_col_name = 'Rechtswert'
y_col_name = 'Hochwert'

# def epsg wgs84 and utm32 for coordinates conversion
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

# threshold for CDF, consider only above thr, below is P0
ppt_thr = .5
max_ppt_thr = 100.

# till 1 day '5min', '10min', '15min', '30min',
aggregation_frequencies = ['60min',
                           '90min', '120min', '180min', '240min',
                           '360min', '480min', '720min', '1440min']
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


def get_for_netatmo_nearest_dwd_station(first_stn_id, dwd_coords_df_file,
                                        netatmo_coords_df_file,
                                        x_col_name, y_col_name,
                                        lon_col_name, lat_col_name):
    ''' Find for one netatmo station, the closest DWD neibhouring station'''
    # read df coordinates and get station ids, and x, y values
    _, x_vals, y_vals, stn_ids = get_dwd_stns_coords(
        dwd_coords_df_file, x_col_name, y_col_name)

    in_netatmo_coords_df, _, _, _ = get_netatmo_stns_coords(
        netatmo_coords_df_file, lon_col_name, lat_col_name)

    first_stn_id = first_stn_id.replace('_', ':')
    lon_stn = in_netatmo_coords_df.loc[first_stn_id, lon_col_name]
    lat_stn = in_netatmo_coords_df.loc[first_stn_id, lat_col_name]

    x_stn, y_stn = convert_coords_fr_wgs84_to_utm32_(wgs82, utm32,
                                                     lon_stn, lat_stn)
    # make a tupples from the coordinates
    coords_tuples = np.array([(x, y) for x, y in zip(x_vals, y_vals)])

    # create a tree from coordinates
    points_tree = spatial.cKDTree(coords_tuples)

    distances, indices = points_tree.query([x_stn, y_stn], k=2)

    coords_nearest_nbr = coords_tuples[indices[1]]
    stn_near = str(stn_ids[indices[1]])
    distance_near = distances[1]

    return coords_nearest_nbr, stn_near, distance_near

#==============================================================================
#
#==============================================================================


def compare_cdf_two_stns(netatmo_ppt_df_file, path_to_ppt_hdf_data):
    HDF52 = HDF5(infile=path_to_ppt_hdf_data)

    in_netatmo_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                     index_col=0, sep=';',
                                     parse_dates=True,
                                     infer_datetime_format=True,
                                     engine='c')
    netatmo_stns_ids = in_netatmo_stns_df.columns

    for stn_id in netatmo_stns_ids[200:]:
        print('First Netatmo Stn Id is', stn_id)

        try:
            idf1 = in_netatmo_stns_df.loc[:, stn_id]
            idf1.dropna(axis=0, inplace=True)
            idf1 = idf1[idf1 < max_ppt_thr]
            _, dwd_stn_near, distance_near = get_for_netatmo_nearest_dwd_station(
                stn_id, coords_dwd_df_file, coords_netatmo_df_file,
                x_col_name, y_col_name, lon_col_name, lat_col_name)

            idf2 = HDF52.get_pandas_dataframe(ids=[dwd_stn_near])
            idf2 = idf2[idf2 < max_ppt_thr]
            print('Second DWD Stn Id is', dwd_stn_near)
            for tem_freq in aggregation_frequencies:
                print('Aggregation is: ', tem_freq)
                df_resample1 = resampleDf(data_frame=idf1,
                                          temp_freq=tem_freq)
                df_resample2 = resampleDf(data_frame=idf2,
                                          temp_freq=tem_freq)
#                 df2 = select_df_within_period(df_resample2,
#                                               df_resample1.index[0],
#                                               df_resample1.index[-1])
#                 df1 = select_df_within_period(df_resample1,
#                                               df2.index[0],
#                                               df2.index[-1])
                idx_common = df_resample1.index.intersection(
                    df_resample2.index)
                df_common1 = df_resample1.loc[idx_common]
                df_common2 = df_resample2.loc[idx_common]
                if (df_common1.values.shape[0] > 0 and
                        df_common2.values.shape[0] > 0):

                    try:
                        plt_bar_plot_2_stns(stn_id, dwd_stn_near, distance_near,
                                            df_common1, df_common2,
                                            tem_freq,
                                            out_save_dir)
                        plt_scatter_plot_2_stns(stn_id, dwd_stn_near, distance_near,
                                                df_common1, df_common2,
                                                tem_freq,
                                                out_save_dir)
                        plot_end_tail_cdf_2_stns(stn_id, dwd_stn_near, distance_near,
                                                 df_common1, df_common2,
                                                 tem_freq, ppt_thr,
                                                 out_save_dir)
                        plot_ranked_stns(stn_id, dwd_stn_near, distance_near,
                                         df_common1, df_common2,
                                         tem_freq,
                                         out_save_dir)
                    except Exception as msg:
                        print('error while plotting', msg, tem_freq)
                        continue
                else:
                    print('empty df')
                    continue
                # break
        except Exception as msg:
            print(msg)

        break


if __name__ == '__main__':
    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    compare_cdf_two_stns(path_to_ppt_netatmo_data,
                         path_to_ppt_hdf_data)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
