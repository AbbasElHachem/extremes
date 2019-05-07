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

from b_get_data import HDF5
# TODO: run script only if df does not exist, MAKE ME FASTER


path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')


out_save_dir = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_08')

if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)

# till 6 hours
aggregation_frequencies = ['10min', '15min', '30min', '60min', '90min',
                           '120min', '180min', '240min', '300min', '360min']
#==============================================================================
# look at different aggregations
#==============================================================================


def resampleDf(data_frame, df_sep_,
               temp_freq,  temp_shift, label_shift,
               out_save_dir=None,
               fillnan=False, df_save_name=None):
    ''' sample DF based on freq and time shift and label shift '''

    df_ = data_frame.copy()
    df_res = df_.resample(temp_freq,
                          label='right',
                          closed='right',
                          loffset=label_shift,
                          base=temp_shift).sum()
    # used to make ppt data same as radolan data UTC
    df_res.index = df_res.index.tz_localize('UTC').tz_convert('Etc/GMT+1')

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


def build_edf_fr_vals(ppt_data, ppt_thr):
    # Construct EDF
    ''' construct empirical distribution function given data values '''
    p0 = calculate_probab_ppt_below_thr(ppt_data, ppt_thr)

    ppt_data_abv_thr = ppt_data[ppt_data > ppt_thr]
    ppt_data_abv_thr = ppt_data
    data_sorted = np.sort(ppt_data_abv_thr, axis=0)[::-1]

    x0 = np.squeeze(data_sorted)[::-1]
    y0 = (np.arange(data_sorted.size) / len(data_sorted))
#     np.append(x0, 0)
    return x0, y0, p0


#==============================================================================
#
#==============================================================================
import matplotlib.pyplot as plt


def compare_cdf_two_stns(stns_ids):
    for iid in stns_ids:
        try:
            idf1 = HDF52.get_pandas_dataframe(ids=[iid])
            idf1.dropna(axis=0, inplace=True)
        except Exception as msg:
            print(msg)
            continue
        x0, y0, p0 = build_edf_fr_vals(idf1.values, 1)
        ids2 = np.array(list(filter(lambda x: x != iid, ids)))
        count_all_stns = len(ids2)
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

        return idf1, x0, y0, p0


if __name__ == '__main__':
    HDF52 = HDF5(infile=path_to_ppt_hdf_data)
    ids = HDF52.get_all_ids()
    compare_cdf_two_stns(ids)
