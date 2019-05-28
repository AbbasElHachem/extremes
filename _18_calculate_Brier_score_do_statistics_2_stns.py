# !/usr/bin/env python.
# -*- coding: utf-8 -*-


"""
Look at pairs of stations (NetAtmo-DWD) 

Select nearest stations, for one DWD station, aggregate to daily
and disaggregate according to relative values of NetAtmo station

Check if this gives reasonable Temporal structure or not 

Save the daily DWD and NetAtmo sums in a dataframe, use it to find 
scaling factor between the values
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
import datetime
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears
from dateutil import parser, rrule


from _00_additional_functions import (resample_intersect_2_dfs, resampleDf,
                                      calculate_probab_ppt_below_thr)


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


distance_matrix_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                           r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD.csv')
assert os.path.exists(distance_matrix_df_file), 'wrong Distance MTX  file'

out_save_dir_orig = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\statistics_Netatmo_Dwd')

if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)


max_ppt_thr = 100.
ppt_thrs = [0.5, 1, 2, 5]
aggregation_frequencies = ['5min', '10min', '15min', '30min', '60min', '90min',
                           '120min', '180min', '240min',  '360min',
                           '480min', '720min', '1440min']

#==============================================================================
#
#==============================================================================


def calculate_brier_score(stn_1_id,  # id of first station
                          stn_2_id,  # id of second station
                          df_1,  # df first station
                          df_2,  # df second station
                          min_ppt_thr  # min ppt thr, above is rain
                          ):
    if not isinstance(df_1, pd.DataFrame):
        df_1 = pd.DataFrame(data=df_1.values,
                            index=df_1.index,
                            columns=[stn_1_id])

    if not isinstance(df_2, pd.DataFrame):
        df_2 = pd.DataFrame(data=df_2.values,
                            index=df_2.index,
                            columns=[stn_2_id])

    df_brier = pd.DataFrame(index=df_1.index,
                            data=np.zeros(shape=(df_1.shape[0], 2)),
                            columns=[stn_1_id, stn_2_id])

    df_brier[stn_1_id] = df_1[df_1[stn_1_id] > min_ppt_thr]
    df_brier[stn_1_id][df_brier[stn_1_id] > min_ppt_thr] = 1
    df_brier[stn_1_id].fillna(0, inplace=True)

    df_brier[stn_2_id] = df_2[df_2[stn_2_id] > min_ppt_thr]
    df_brier[stn_2_id][df_brier[stn_2_id] > min_ppt_thr] = 1
    df_brier[stn_2_id].fillna(0, inplace=True)

    df_brier['Score_diff'] = np.square(df_brier[stn_1_id] - df_brier[stn_2_id])

    brier_score = df_brier['Score_diff'].sum() / df_brier.shape[0]
    return brier_score
#==============================================================================
#
#==============================================================================


def calculate_brier_score_2_stns_per_temp_freq(netatmo_ppt_df_file,
                                               path_to_ppt_hdf_data,
                                               distance_matrix_df_file,
                                               aggregation_frequencies_list,
                                               ppt_thrs_list,
                                               out_dir):
    HDF52 = HDF5(infile=path_to_ppt_hdf_data)

    in_netatmo_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                     index_col=0, sep=';',
                                     parse_dates=True,
                                     infer_datetime_format=True,
                                     engine='c')
    netatmo_stns_ids = in_netatmo_stns_df.columns

    in_df_distance_netatmo_dwd = pd.read_csv(distance_matrix_df_file,
                                             sep=';', index_col=0)

    for stn_id in netatmo_stns_ids:

        print('First Netatmo Stn Id is', stn_id)

        try:
            idf1 = in_netatmo_stns_df.loc[:, stn_id]
            idf1.dropna(axis=0, inplace=True)
            idf1 = idf1[idf1 < max_ppt_thr]

            distances_to_stn1 = in_df_distance_netatmo_dwd.loc[stn_id, :]
            sorted_distances = distances_to_stn1.sort_values(
                ascending=True)

            min_dist = sorted_distances.values[0]
            if min_dist <= 5000:

                stn_2_id = sorted_distances.index[0]

                idf2 = HDF52.get_pandas_dataframe(ids=[stn_2_id])
                idf2 = idf2[idf2 < max_ppt_thr]

                print('Second DWD Stn Id is', stn_2_id)

                for ppt_thr in ppt_thrs_list:
                    df_statistics_with_agg = pd.DataFrame(
                        index=aggregation_frequencies_list)

                    for temp_freq in aggregation_frequencies_list:
                        print('Temp frequency is', temp_freq)
                        df_netatmo, df_dwd = resample_intersect_2_dfs(
                            idf1,
                            idf2,
                            temp_freq)

                        if (df_netatmo.values.shape[0] > 1000):
                            raise Exception
                            print('getting statistics')
                            brier_score = calculate_brier_score(stn_1_id=stn_id,
                                                                stn_2_id=stn_2_id,
                                                                df_1=df_netatmo,
                                                                df_2=df_dwd,
                                                                min_ppt_thr=ppt_thr)
                            # calculate correlations (pearson and spearman)
                            corr = pears(df_dwd.values.ravel(),
                                         df_netatmo.values.ravel())[0]
                            rho = spr(df_dwd.values.ravel(),
                                      df_netatmo.values.ravel())[0]
                            p01 = calculate_probab_ppt_below_thr(
                                df_netatmo.values, ppt_thr)
                            p02 = calculate_probab_ppt_below_thr(
                                df_dwd.values, ppt_thr)

                            df_statistics_with_agg.loc[
                                temp_freq,
                                'Brier Score'] = np.round(brier_score, 3)

                            df_statistics_with_agg.loc[
                                temp_freq,
                                'Pearson Correlation'] = np.round(corr, 3)

                            df_statistics_with_agg.loc[
                                temp_freq,
                                'Spearman Correlation'] = np.round(rho, 2)

                            df_statistics_with_agg.loc[
                                temp_freq,
                                'p0 Netatmo'] = np.round(p01, 2)

                            df_statistics_with_agg.loc[
                                temp_freq,
                                'p0 DWD'] = np.round(p02, 2)
                        else:
                            print('not enough data')

#                     print(df_statistics_with_agg.values.shape[0])

#                     if df_statistics_with_agg.values.shape[0] > 0:
#                         df_statistics_with_agg.to_csv(
#                             os.path.join(out_dir,
#                                          'statistics_netatmo_%s_dwd_%s_ppt_thr_%0.1f_.csv'
#                                          % (stn_id, stn_2_id, ppt_thr)), sep=';',
#                             float_format='0.3f')

        except Exception as msg:
            print(msg)

#         break


if __name__ == '__main__':
    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    calculate_brier_score_2_stns_per_temp_freq(path_to_ppt_netatmo_data,
                                               path_to_ppt_hdf_data,
                                               distance_matrix_df_file,
                                               aggregation_frequencies,
                                               ppt_thrs,
                                               out_save_dir_orig)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
