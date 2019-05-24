# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""Plot for every extreme event at a station the neighbouring stations

Parameters
----------
Path to where simulatenous events dataframes are present

Returns
---------

Plot for every shifted time (+-5min) what was happening at the 
center station and surrounding stations (radius = 30Km)
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================
import os
import timeit
import time


import scipy.spatial as spatial

from adjustText import adjust_text

from matplotlib import rc
from matplotlib import rcParams

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapefile
plt.ioff()
# get_ipython().run_line_magic('matplotlib', 'inline')

rc('font', size=16)
rc('font', family='serif')
rc('axes', labelsize=20)
rcParams['axes.labelpad'] = 35

from _00_additional_functions import (convert_coords_fr_wgs84_to_utm32_,
                                      list_all_full_path,
                                      pltcolor)

from _02_get_data_simultaneous_stns import get_events_stn_data


# plt.style.use('fast')
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'axes.labelsize': 12})

#==============================================================================
#
#==============================================================================
path_to_dfs_simultaneous_events = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\thr8'

out_plots_dir = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\thr8Plots_DWD'

path_to_ppt_coords_data = (r'X:\exchange\ElHachem'
                           r'\niederschlag_deutschland'
                           r'\1993_2016_5min_merge_nan.csv')

path_to_shpfile = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                   r'\shp_file_germany\DEU_adm1.shp')
# in station df, to get coordinates
xcoord_name = 'Rechtswert'
ycoord_name = 'Hochwert'

# def epsg wgs84 and utm32 for coordinates conversion
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

# used when plotting to make different markers per time and values
markers = ['o', '.', ',', '2', '+', 'v', '^', '<', '>', 's', 'd', '1', '8',
           '1', 'd', 's', '>', '<', '^', 'v', '+',  '2', ',', '.', 'o']

markers_time_dict = {i: m for i, m in zip(np.arange(-60, 61, 5), markers)}
#==============================================================================
#
#==============================================================================


def calculate_distance_between_two_positions(x0, y0, x1, y1):
    ''' calculate angle between two successive positions in Km'''
    distance = np.round(
        np.sqrt(np.square(x1 - x0) + np.square(y1 - y0)) / 1000, 2)
    return distance

#==============================================================================
#
#==============================================================================


def direction_lookup(destination_x, origin_x, destination_y, origin_y):
    ''' calculate angle and direction between two stations'''
    deltaX = destination_x - origin_x
    deltaY = destination_y - origin_y
    degrees_temp = np.math.atan2(deltaX, deltaY) / np.math.pi * 180
    if degrees_temp < 0:
        degrees_final = 360 + degrees_temp
    else:
        degrees_final = degrees_temp
    compass_brackets = ["N", "NE", "E", "SE", "S", "SW", "W", "NW", "N"]
    compass_lookup = round(degrees_final / 45)
    return compass_brackets[compass_lookup], np.round(degrees_final, 3)

#==============================================================================
#
#==============================================================================


def plot_stations_in_convex_hull(path_to_events,
                                 path_stns_coords,
                                 xcoords_name,
                                 ycoords_name,
                                 shp_de,
                                 radius=30000):
    '''fct to plot station event data and all other stations within +-60min'''
    colors_dict = pltcolor()
    # get all events, extracted from script _00_
    dfs_data_lst = list_all_full_path('.csv', path_to_events)

    shp_de = shapefile.Reader(shp_de)

    for evt_idx, event in enumerate(dfs_data_lst):

        (stn_one_id, ppt_stn_one, stns_2_ids,
         stn_one_xcoords, stn_one_ycoords,
         stns_2_xcoords, stns_2_ycoords,
         event_date, stns_2_ids_vals_dict,
         ppt_thr) = get_events_stn_data(event,
                                        path_stns_coords,
                                        xcoords_name,
                                        ycoords_name)

        out_path = os.path.join(out_plots_dir, str(int(stn_one_id)))
        if not os.path.exists(out_path):
            os.mkdir(out_path)

        out_event_path = os.path.join(out_path, str(int(evt_idx)))
        if not os.path.exists(out_event_path):
            os.mkdir(out_event_path)

        coords_tuples = np.array([(x2, y2) for x2, y2
                                  in zip(stns_2_xcoords, stns_2_ycoords)])

        points_tree = spatial.cKDTree(coords_tuples)

        # This finds the index of all points within distance 1 of [1.5,2.5].
        idxs_neighbours = points_tree.query_ball_point(
            np.array((stn_one_xcoords, stn_one_ycoords)), radius)

        stns2_ids_in_circle = np.empty(shape=len(idxs_neighbours))
        stns2_xcoords_in_circle = np.empty(shape=len(idxs_neighbours))
        stns2_ycoords_in_circle = np.empty(shape=len(idxs_neighbours))

        for i, ix_nbr in enumerate(idxs_neighbours):
            stns2_ids_in_circle[i] = stns_2_ids[ix_nbr]
            stns2_xcoords_in_circle[i] = stns_2_xcoords[ix_nbr]
            stns2_ycoords_in_circle[i] = stns_2_ycoords[ix_nbr]

        stns_2_ids_vals_dict_in_circle = {}
        for k in stns2_ids_in_circle:
            if k in stns_2_ids_vals_dict.keys():
                stns_2_ids_vals_dict_in_circle[k] = stns_2_ids_vals_dict[k]

        time_vals = np.arange(-60, 61, 5)
        df_out = pd.DataFrame(
            index=time_vals, columns=stns_2_ids_vals_dict_in_circle.keys())
        save_event_time = event_date.replace(
            ':', '_').replace(' ', '_').replace('-', '_')
        for timeval in time_vals:
            print(timeval)
            # plot all other simultaneous stations
            texts = []

            fig = plt.figure(figsize=(20, 20), dpi=75)
            ax = fig.add_subplot(111)
            for shape_ in shp_de.shapeRecords():
                lon = [i[0]
                       for i in shape_.shape.points[:][::-1]]
                lat = [i[1]
                       for i in shape_.shape.points[:][::-1]]

                x0, y0 = convert_coords_fr_wgs84_to_utm32_(
                    wgs82, utm32, lon, lat)

                ax.scatter(x0, y0, marker='.',
                           c='lightgrey', alpha=0.5, s=2)

                # plot first station
                ax.scatter(stn_one_xcoords,
                           stn_one_ycoords,
                           c='red',
                           marker='X',
                           s=50,
                           label=('Stn %s Ppt %0.2f mm'
                                  % (str(int(stn_one_id)),
                                     ppt_stn_one)))

            for i, stn2_id in enumerate(stns2_ids_in_circle):
                stn2_xcoord = stns2_xcoords_in_circle[i]
                stn2_ycoord = stns2_ycoords_in_circle[i]
                time_idx_val = stns_2_ids_vals_dict_in_circle[stn2_id][timeval]
                # get distance and angle
                dist = calculate_distance_between_two_positions(
                    stn_one_xcoords, stn_one_ycoords, stn2_xcoord, stn2_ycoord)
                orient, _ = direction_lookup(
                    stn2_xcoord, stn_one_xcoords, stn2_ycoord, stn_one_ycoords)

                df_out.loc['Distanz (Km)', stn2_id] = dist
                df_out.loc['Richtung', stn2_id] = orient
                if len(time_idx_val) > 0:
                    df_out.loc[timeval, stn2_id] = np.float(time_idx_val[0])
#
                    ax.scatter(stn2_xcoord,
                               stn2_ycoord,
                               c=colors_dict[timeval][0],
                               marker=markers_time_dict[timeval],
                               s=45,
                               label=('Stn %s Ppt %0.2f at %0.0f min'
                                      % (str(int(stn2_id)),
                                         time_idx_val[0],
                                         timeval)))

#                     val_float_format = '% 0.1f mm %0.0f min' % (
#                         time_idx_val[0], timeval)
                    val_float_format = '% 0.1f mm' % (time_idx_val[0])

                    texts.append(ax.text(stn2_xcoord,
                                         stn2_ycoord,
                                         val_float_format,
                                         color='k'))
                    # color=colrs_dict[timeval][0]))

            ppt_stn_one_float_format = '% 0.2f' % ppt_stn_one
            texts.append(ax.text(stn_one_xcoords,
                                 stn_one_ycoords,
                                 ppt_stn_one_float_format))

            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle='->', color='red', lw=0.25))

            ax.set_title('Event_at_station_%s_ppt_thr_%smm_at_%s_at_time_%dmin'
                         % (str(int(stn_one_id)), ppt_thr,
                             event_date, timeval))

            ax.grid(alpha=0.25)

            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

            ax.set_xlim([250000, 950000])
            ax.set_ylim([5200000, 6110000])
            ax.set_aspect(1.0)

            plt.tight_layout()

            plt.savefig(
                os.path.join(out_event_path,
                             'station_%s_ppt_thr_%smm_at_%s_2_%dmin_.png'
                             % (str(int(stn_one_id)),
                                ppt_thr, save_event_time,
                                timeval)),
                frameon=True, papertype='a4',
                bbox_inches='tight', pad_inches=.2)
            plt.close()
            print('Finished plotting')
        df_out.to_csv(
            os.path.join(out_event_path,
                         'station_%s_ppt_thr_%smm_at_%s_2_%dmin_.csv'
                         % (str(int(stn_one_id)),
                            ppt_thr, save_event_time,
                            timeval)), sep=';', float_format='%.2f')

    return df_out
#==============================================================================
#
#==============================================================================


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    df_out = plot_stations_in_convex_hull(path_to_dfs_simultaneous_events,
                                          path_to_ppt_coords_data, xcoord_name,
                                          ycoord_name, path_to_shpfile,
                                          30000)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
