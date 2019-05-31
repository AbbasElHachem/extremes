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


from _00_additional_functions import (resample_intersect_2_dfs)

from _09_aggregate_plot_compare_2_DWD_stns import (plt_bar_plot_2_stns,
                                                   plt_scatter_plot_2_stns,
                                                   plot_end_tail_cdf_2_stns,
                                                   plot_normalized_ranked_stns,
                                                   plot_normalized_sorted_ranked_stns,
                                                   plot_sorted_stns_vals,
                                                   plot_p0_as_a_sequence_two_stns,
                                                   plot_contingency_tables_as_a_sequence_two_stns)
from b_get_data import HDF5
#==============================================================================
#
#==============================================================================

path_to_ppt_hourly_netatmo_data = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(
    path_to_ppt_hourly_netatmo_data), 'wrong NETATMO Ppt file'

path_to_ppt_5min_netatmo_data = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_5min_stns_combined_.csv')
assert os.path.exists(
    path_to_ppt_5min_netatmo_data), 'wrong NETATMO Ppt file'

path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')
assert os.path.exists(path_to_ppt_hdf_data), 'wrong NETATMO Ppt file'


distance_matrix_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                           r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD.csv')
assert os.path.exists(distance_matrix_df_file), 'wrong Distance MTX  file'

# out_save_dir_orig = (
#     r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_DWD_NetAtmo_hourly')

out_save_dir_orig = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\cdf_plots_DWD_NetAtmo_5min')

if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)


# threshold for CDF, consider only above thr, below is P0
ppt_thr = .5
ppt_thr2 = 0.5


ppt_thrs_list = [0.5, 1., 2]

max_ppt_thr = 100.

# till 1 day '5min', '10min', '15min', '30min',
aggregation_frequencies = ['5min', '10min', '15min', '30min', '60min',
                           '90min', '120min', '180min', '240min',
                           '360min', '480min', '720min', '1440min']
#==============================================================================
#
#==============================================================================


def compare_cdf_two_stns(netatmo_ppt_df_file, path_to_ppt_hdf_data,
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
            assert stn_id in in_df_distance_netatmo_dwd.index, 'check dist mtx'
            distances_to_stn1 = in_df_distance_netatmo_dwd.loc[stn_id, :]
            sorted_distances = distances_to_stn1.sort_values(
                ascending=True)

            min_dist = sorted_distances.values[0]
            if min_dist <= 5000:

                stn_2_id = sorted_distances.index[0]

                idf2 = HDF52.get_pandas_dataframe(ids=[stn_2_id])
                idf2 = idf2[idf2 < max_ppt_thr]
                print('Second DWD Stn Id is', stn_2_id,
                      'distance is ', min_dist)
                if (idf1.values.shape[0] > 1000) and (idf2.values.shape[0] > 1000):
                    out_save_dir = os.path.join(out_save_dir_orig,
                                                '%s_%s' % (stn_id, stn_2_id))
                    if not os.path.exists(out_save_dir):
                        os.mkdir(out_save_dir)

                    for tem_freq in aggregation_frequencies:
                        print('Aggregation is: ', tem_freq)

                        df_common1, df_common2 = resample_intersect_2_dfs(idf1,
                                                                          idf2,
                                                                          tem_freq)

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

                                plot_normalized_sorted_ranked_stns(stn_id,
                                                                   stn_2_id,
                                                                   min_dist,
                                                                   df_common1,
                                                                   df_common2,
                                                                   tem_freq,
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
                            break
                    plot_p0_as_a_sequence_two_stns(stn_id,
                                                   stn_2_id,
                                                   min_dist,
                                                   ppt_thrs_list,
                                                   idf1,
                                                   idf2,
                                                   aggregation_frequencies,
                                                   out_save_dir)

                    plot_contingency_tables_as_a_sequence_two_stns(stn_2_id,
                                                                   stn_id,
                                                                   min_dist,
                                                                   ppt_thrs_list,
                                                                   idf2,
                                                                   idf1,
                                                                   aggregation_frequencies,
                                                                   out_save_dir)
                else:
                    print('Station is near but dont have enough data')
            else:
                print('station is not near looking for another station')
                continue

        except Exception as msg:
            print(msg)
#         break


if __name__ == '__main__':
    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    compare_cdf_two_stns(path_to_ppt_5min_netatmo_data,
                         path_to_ppt_hdf_data,
                         distance_matrix_df_file)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
