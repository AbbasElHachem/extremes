# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
For every station, and every extreme event plot the station and 
all other stations where simulatneously (+- 60min) an extreme 
event was occurring
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

import os
import timeit
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy import spatial


from _00_additional_functions import (list_all_full_path, resampleDf)
from _06_get_data_simultaneous_Netatmo_events import get_netatmo_events_stn_data

from _09_aggregate_do_cdf_compare_2_DWD_stns import (plt_bar_plot_2_stns,
                                                     plt_scatter_plot_2_stns,
                                                     plot_end_tail_cdf_2_stns,
                                                     plot_ranked_stns)

from _10_aggregate_do_cdf_compare_2_NetAtmo_stns import aggregation_frequencies

# plt.style.use('fast')
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 12})

#==============================================================================
#
#==============================================================================
path_to_dfs_simultaneous_events = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\thr24NetAtmo'

path_to_netatmo_coords_df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
                                  r"\rain_bw_1hour"
                                  r"\netatmo_bw_1hour_coords.csv")

path_to_netatmo_data_df_all = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
                               r"\ppt_all_netatmo_hourly_stns_combined_.csv")

out_plots_dir = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmoPlotsExtremes"

if not os.path.exists(out_plots_dir):
    os.mkdir(out_plots_dir
             )
# in station df, to get coordinates
xcoord_name = ' lon'
ycoord_name = ' lat'


max_ppt_thr = 60.
p0_ppt_thr = 1.
#==============================================================================
#
#==============================================================================


def plot_netatmot_coordinates(path_to_events,
                              path_stns_coords,
                              path_to_netatmo_ppt_df,
                              xcoords_name,
                              ycoords_name,
                              radius,
                              out_save_dir):
    '''fct to plot station event data and all other stations within +-60min'''

    in_netatmo_stns_df = pd.read_csv(path_to_netatmo_ppt_df,
                                     index_col=0, sep=';',
                                     parse_dates=True,
                                     infer_datetime_format=True,
                                     engine='c')

    # get all events, extracted from script _00_
    dfs_data_lst = list_all_full_path('.csv', path_to_events)

    for i, event in enumerate(dfs_data_lst[2:]):

        try:
            (stn_one_id, ppt_stn_one, stns_2_ids,
             stn_one_xcoords, stn_one_ycoords,
             stns_2_xcoords, stns_2_ycoords,
             event_date, stns_2_ids_vals_dict,
             ppt_thr) = get_netatmo_events_stn_data(event,
                                                    path_stns_coords,
                                                    xcoords_name,
                                                    ycoords_name)

            save_event_time = event_date.replace(
                ':', '_').replace(' ', '_').replace('-', '_')

            if not os.path.exists(
                os.path.join(out_save_dir,
                             'station_%s_ppt_thr_%smm_at_%s_2.png'
                             % (str(stn_one_id),
                                ppt_thr, save_event_time))):

                coords_tuples = np.array([(x2, y2) for x2, y2
                                          in zip(stns_2_xcoords,
                                                 stns_2_ycoords)])

                points_tree = spatial.cKDTree(coords_tuples)

                # This finds the index of all points within distance radius
                idxs_neighbours = points_tree.query_ball_point(
                    np.array((stn_one_xcoords, stn_one_ycoords)), radius)

                stns2_ids_in_circle = []
                stns2_xcoords_in_circle = np.empty(shape=len(idxs_neighbours))
                stns2_ycoords_in_circle = np.empty(shape=len(idxs_neighbours))

                for i, ix_nbr in enumerate(idxs_neighbours):
                    stns2_ids_in_circle.append(stns_2_ids[ix_nbr])
                    stns2_xcoords_in_circle[i] = stns_2_xcoords[ix_nbr]
                    stns2_ycoords_in_circle[i] = stns_2_ycoords[ix_nbr]

                stns_2_ids_vals_dict_in_circle = {}
                for k in stns2_ids_in_circle:
                    if k in stns_2_ids_vals_dict.keys():
                        stns_2_ids_vals_dict_in_circle[k] = stns_2_ids_vals_dict[k]

                print('First Stn Id is', stn_one_id)
                try:
                    idf1 = in_netatmo_stns_df.loc[:, stn_one_id]
                    idf1.dropna(axis=0, inplace=True)
                    idf1 = idf1[idf1 < max_ppt_thr]

                    for i, stn2_id in enumerate(stns2_ids_in_circle):

                        idf2 = in_netatmo_stns_df.loc[:, stn2_id]

                        idf2 = idf2[idf2 < max_ppt_thr]
                        print('Second Stn Id is', stn2_id)
                        for tem_freq in aggregation_frequencies:
                            print('Aggregation is: ', tem_freq)

                            df_resample1 = resampleDf(data_frame=idf1,
                                                      temp_freq=tem_freq)
                            df_resample2 = resampleDf(data_frame=idf2,
                                                      temp_freq=tem_freq)
#                             df_resample1 = df_resample1[df_resample1 >= 0.0]

                            idx_common = df_resample1.index.intersection(
                                df_resample2.index)
                            df_common1 = df_resample1.loc[idx_common]
                            df_common2 = df_resample2.loc[idx_common]
                            if (df_common1.values.shape[0] > 0 and
                                    df_common2.values.shape[0] > 0):
                                distance_near = radius
                                try:
                                    plt_bar_plot_2_stns(stn_one_id,
                                                        stn2_id,
                                                        distance_near,
                                                        df_common1,
                                                        df_common2,
                                                        tem_freq,
                                                        out_save_dir)
                                    plt_scatter_plot_2_stns(stn_one_id,
                                                            stn2_id,
                                                            distance_near,
                                                            df_common1,
                                                            df_common2,
                                                            tem_freq,
                                                            out_save_dir)
                                    plot_end_tail_cdf_2_stns(stn_one_id,
                                                             stn2_id,
                                                             distance_near,
                                                             df_common1,
                                                             df_common2,
                                                             tem_freq,
                                                             p0_ppt_thr,
                                                             out_save_dir)
                                    plot_ranked_stns(stn_one_id, stn2_id,
                                                     distance_near,
                                                     df_common1, df_common2,
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
                #==============================================================

            else:
                print('plots already exists')
                continue
        except Exception as msg:
            print(msg)
            continue
    return


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    plot_netatmot_coordinates(path_to_events=path_to_dfs_simultaneous_events,
                              path_stns_coords=path_to_netatmo_coords_df_file,
                              path_to_netatmo_ppt_df=path_to_netatmo_data_df_all,
                              xcoords_name=xcoord_name,
                              ycoords_name=ycoord_name,
                              radius=30000,
                              out_save_dir=out_plots_dir)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
