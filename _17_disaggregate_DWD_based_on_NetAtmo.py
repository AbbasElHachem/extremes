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
import matplotlib.pyplot as plt

import datetime
from _00_additional_functions import (resample_intersect_2_dfs, resampleDf)

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
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\disaggregate_DWD_NetAtmo')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)


max_ppt_thr = 100.

#==============================================================================
#
#==============================================================================


def construct_netatmo_dwd_daily_dfs(netatmo_ppt_df_file,
                                    path_to_ppt_hdf_data,
                                    distance_matrix_df_file):
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
                raise Exception
                stn_2_id = sorted_distances.index[0]

                idf2 = HDF52.get_pandas_dataframe(ids=[stn_2_id])
                idf2 = idf2[idf2 < max_ppt_thr]

                print('Second DWD Stn Id is', stn_2_id)

                df_combined = pd.DataFrame(columns=[stn_id, stn_2_id])

                df_netatmo_hourly, df_dwd_hourly = resample_intersect_2_dfs(
                    idf1,
                    idf2,
                    '60min')
                if (df_netatmo_hourly.values.shape[0] > 1000):
                    raise Exception
                    df_dwd_daily = resampleDf(df_dwd_hourly, '1440min')
                    df_netatmo_daily = resampleDf(df_netatmo_hourly, '1440min')

                    df_disaggregated = pd.DataFrame(data=np.empty(shape=(df_dwd_hourly.shape[0])),
                                                    index=df_dwd_hourly.index)

                    for idx, val in zip(df_dwd_hourly.index,
                                        df_dwd_hourly.values):
                        daily_idx = datetime.date(idx.year, idx.month, idx.day)
                        if val == 0:
                            print('Val is 0')
                            df_disaggregated.loc[idx] = val

                        if val != 0:
                            if daily_idx in df_dwd_daily.index:
                                print('val is', val, 'at ', daily_idx)
                                sum_daily_dwd = df_dwd_daily.loc[daily_idx].values
                                sum_daily_netatmo = df_netatmo_daily.loc[daily_idx]
                                # TODO FIX ME
                                print('sum dwd is ', sum_daily_dwd)
                                print('sum netatmo is', sum_daily_netatmo)
                                if sum_daily_dwd == 0:
                                    ppt_disagg = 0
                                else:
                                    ppt_disagg = (
                                        val * sum_daily_netatmo) / sum_daily_dwd

                                print('disaggregated ppt is', ppt_disagg)
                                df_disaggregated.loc[idx] = ppt_disagg
                            else:
                                df_disaggregated.loc[idx] = np.nan
                            # break
                    plt.ioff()
                    plt.plot(df_disaggregated.index,
                             df_disaggregated.values, c='r', alpha=0.25)
                    plt.plot(df_dwd_hourly.index,
                             df_dwd_hourly.values, c='b', alpha=0.25)
                    plt.show()
                else:
                    print('empty df')
                    continue

        except Exception as msg:
            print(msg)
#         break


if __name__ == '__main__':
    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    construct_netatmo_dwd_daily_dfs(path_to_ppt_netatmo_data,
                                    path_to_ppt_hdf_data,
                                    distance_matrix_df_file)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
