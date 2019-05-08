# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Look at pairs of stations (either DWD-DWD, DWD-Netatmo, Netatmo-Netatmo) 

For different aggregations (15min, 30min, 60min, ..., 6hours)
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

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy.spatial as spatial
from scipy import stats

from b_get_data import HDF5


path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')

coords_df_file = r'X:\exchange\ElHachem\niederschlag_deutschland\1993_2016_5min_merge_nan.csv'
assert os.path.exists(coords_df_file), 'wrong DWD coords file'

out_save_dir = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_08')


if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)

x_col_name = 'Rechtswert'
y_col_name = 'Hochwert'

# till 6 hours
aggregation_frequencies = ['10min', '15min', '30min', '60min', '90min',
                           '120min', '180min', '240min', '300min', '360min']
#==============================================================================
# look at different aggregations
#==============================================================================


def resampleDf(data_frame,
               temp_freq,
               temp_shift=0,
               label_shift=None,
               df_sep_=None,
               out_save_dir=None,
               fillnan=False,
               df_save_name=None):
    ''' sample DF based on freq and time shift and label shift '''

    df_ = data_frame.copy()
    df_res = df_.resample(temp_freq,
                          label='right',
                          closed='right',
                          loffset=label_shift,
                          base=temp_shift).sum()

    if fillnan:
        df_res.fillna(value=0, inplace=True)
    if df_save_name is not None and out_save_dir is not None:
        df_res.to_csv(os.path.join(out_save_dir, df_save_name),
                      sep=df_sep_)
    return df_res

#==============================================================================
#
#==============================================================================


def calculate_probab_ppt_below_thr(ppt_data, ppt_thr):
    ''' calculate probability of values being below threshold '''
    origin_count = ppt_data.shape[0]
    count_below_thr = ppt_data[ppt_data <= ppt_thr].shape[0]
    p0 = np.divide(count_below_thr, origin_count)
    return p0


#==============================================================================
#
#==============================================================================


def build_edf_fr_vals(ppt_data):
    # Construct EDF
    ''' construct empirical distribution function given data values '''
    data_sorted = np.sort(ppt_data, axis=0)[::-1]
    x0 = np.squeeze(data_sorted)[::-1]
    y0 = (np.arange(data_sorted.size) / len(data_sorted))
    return x0, y0


#==============================================================================
#
#==============================================================================


def get_cdf_part_abv_thr(ppt_data, ppt_thr):
    ''' select part of the CDF that is abv ppt thr '''

    p0 = calculate_probab_ppt_below_thr(ppt_data, ppt_thr)

    x0, y0 = build_edf_fr_vals(ppt_data)
    x_abv_thr = x0[x0 > ppt_thr]
    y_abv_thr = y0[np.where(x0 > ppt_thr)]
    assert y_abv_thr[0] == p0, 'something is wrong with probability cal'

    return x_abv_thr, y_abv_thr


def get_dwd_stns_coords(coords_df_file, x_col_name, y_col_name):
    in_coords_df = pd.read_csv(coords_df_file, sep=';',
                               index_col=3, engine='c')
    stn_ids = in_coords_df.index
    x_vals = in_coords_df[x_col_name].values.ravel()
    y_vals = in_coords_df[y_col_name].values.ravel()
    return in_coords_df, x_vals, y_vals, stn_ids


def get_nearest_station(first_stn_id, coords_df_file,
                        x_col_name, y_col_name):
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
    stn_near = stn_ids[indices[1]]
    distance_near = distances[1]

    return coords_nearest_nbr, stn_near, distance_near


def plt_scatter_plot_2_stns(df1, df2):
    plt.scatter(df1.values, df2.values, c='b', alpha=0.5)
    return


def compare_cdf_two_stns(stns_ids):
    for iid in stns_ids:
        try:
            idf1 = HDF52.get_pandas_dataframe(ids=[iid])
            idf1.dropna(axis=0, inplace=True)
            df_resample1 = resampleDf(data_frame=idf1,
                                      temp_freq=aggregation_frequencies[0])

            coords_near, stn_near, distance_near = get_nearest_station(
                iid, coords_df_file, x_col_name, y_col_name)
            assert iid != str(stn_near), 'wrong neighbour selected'
            idf2 = HDF52.get_pandas_dataframe(ids=str(stn_near))

            df_resample2 = resampleDf(data_frame=idf2,
                                      temp_freq=aggregation_frequencies[0])

            idx_common = df_resample1.index.intersection(df_resample2.index)
            # TODO fix plotting
            df_common1 = df_resample1.loc[idx_common, :]  # .rank()
            df_common2 = df_resample2.loc[idx_common, :].rank()

            rank_data = stats.rankdata(df_common1.values)
            plt.scatter(df_common1.values, df_common2.values, c='b', alpha=0.5)
        except Exception as msg:
            print(msg)
            continue

        break
#         x_extremes, y_extremes = get_cdf_part_abv_thr(idf1.values, 1)
#         ids2 = np.array(list(filter(lambda x: x != iid, ids)))
#         count_all_stns = len(ids2)
#         for ii2, iid2 in enumerate(ids2):
#
#             print('Second Station ID is: ', iid2,
#                   ' index is ', ii2,
#                   ' Count of Stns is :', count_all_stns)
#             # read second station
#             try:
#                 idf2 = HDF52.get_pandas_dataframe(ids=[iid2])
#             except Exception as msg:
#                 print(msg)
#                 continue

#         return idf1, x_extremes, y_extremes


if __name__ == '__main__':
    HDF52 = HDF5(infile=path_to_ppt_hdf_data)
    ids = HDF52.get_all_ids()
    compare_cdf_two_stns(ids)
#     coords_nearest_nbr, stn_near = get_nearest_station(coords_df_file)
