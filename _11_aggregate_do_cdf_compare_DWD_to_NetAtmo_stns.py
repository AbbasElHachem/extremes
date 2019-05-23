# !/usr/bin/env python.
# -*- coding: utf-8 -*-


"""
Look at pairs of stations (NetAtmo-DWD) 

For different aggregations (15min, 30min, 60min, ..., 6hours, ..., 1 day)
Calculate P0 (probability that rainfall < 1mm)
Construct Cdf for extremes ( > 1mm) and compare
Calculate the correlation between the ranks on all scales
Find if the similarities or differences in the stations are due to 
a systematic or a random error.
If a systematic error is present (bias) find the correction factor 

Look at pairs of stations (DWD-NetAtmo) 

Construct Contingency Tables for closest pair of stations 
        
             below thr above thr   
below thr    %         %
above thr    %         %

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


from _00_additional_functions import (resampleDf)

from _09_aggregate_do_cdf_compare_2_DWD_stns import (plt_bar_plot_2_stns,
                                                     plt_scatter_plot_2_stns,
                                                     plot_end_tail_cdf_2_stns,
                                                     plot_normalized_ranked_stns,
                                                     plot_sorted_stns_vals)
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

out_save_dir = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_DWD_NetAtmo')


if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)


# threshold for CDF, consider only above thr, below is P0
ppt_thr = .5
ppt_thr2 = 0

max_ppt_thr = 100.

# till 1 day '5min', '10min', '15min', '30min',
aggregation_frequencies = ['60min',
                           '90min', '120min', '180min', '240min',
                           '360min', '480min', '720min', '1440min']
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

    in_df_distance_netatmo_dwd = pd.read_csv(distance_matrix_df_file,
                                             sep=';', index_col=0)

    for stn_id in netatmo_stns_ids[200:]:
        print('First Netatmo Stn Id is', stn_id)

        try:
            idf1 = in_netatmo_stns_df.loc[:, stn_id]
            idf1.dropna(axis=0, inplace=True)
            idf1 = idf1[idf1 < max_ppt_thr]

            distances_to_stn1 = in_df_distance_netatmo_dwd.loc[stn_id, :]
            sorted_distances = distances_to_stn1.sort_values(ascending=True)

            min_dist = sorted_distances.values[0]
            stn_2_id = sorted_distances.index[0]

            idf2 = HDF52.get_pandas_dataframe(ids=[stn_2_id])
            idf2 = idf2[idf2 < max_ppt_thr]
            print('Second DWD Stn Id is', stn_2_id)

            for tem_freq in aggregation_frequencies:
                print('Aggregation is: ', tem_freq)
                df_resample1 = resampleDf(data_frame=idf1,
                                          temp_freq=tem_freq)
                df_resample2 = resampleDf(data_frame=idf2,
                                          temp_freq=tem_freq)

                idx_common = df_resample1.index.intersection(
                    df_resample2.index)
                df_common1 = df_resample1.loc[idx_common]
                df_common2 = df_resample2.loc[idx_common]
                if (df_common1.values.shape[0] > 0 and
                        df_common2.values.shape[0] > 0):
                    #
                    try:
                        plt_bar_plot_2_stns(stn_id,
                                            stn_2_id,
                                            min_dist,
                                            df_common1,
                                            df_common2,
                                            tem_freq,
                                            out_save_dir)

                        plt_scatter_plot_2_stns(stn_id,
                                                stn_2_id,
                                                min_dist,
                                                df_common1,
                                                df_common2,
                                                ppt_thr2,
                                                tem_freq,
                                                out_save_dir)

                        plot_end_tail_cdf_2_stns(stn_id,
                                                 stn_2_id,
                                                 min_dist,
                                                 df_common1,
                                                 df_common2,
                                                 tem_freq,
                                                 ppt_thr,
                                                 out_save_dir)

                        plot_normalized_ranked_stns(stn_id,
                                                    stn_2_id,
                                                    min_dist,
                                                    df_common1,
                                                    df_common2,
                                                    tem_freq,
                                                    out_save_dir)

                        plot_sorted_stns_vals(stn_id,
                                              stn_2_id,
                                              min_dist,
                                              df_common1,
                                              df_common2,
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


if __name__ == '__main__':
    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    compare_cdf_two_stns(path_to_ppt_netatmo_data,
                         path_to_ppt_hdf_data)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
