# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""Gets and prints the spreadsheet's header columns

Parameters
----------
file_loc : str
    The file location of the spreadsheet
print_cols : bool, optional
    A flag used to print the columns to the console (default is False)

Returns
-------
list
    a list of strings representing the header columns
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================
import os
import timeit
import time


import scipy.spatial as spatial

from adjustText import adjust_text
from matplotlib.ticker import LinearLocator
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

import matplotlib as mpl
import matplotlib.colors as mcolors
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

# Function to map the colors as a list from the input list of x variables


def plot_stations_in_convex_hull(path_to_events,
                                 path_stns_coords,
                                 xcoords_name,
                                 ycoords_name,
                                 shp_de,
                                 radius=30000):
    '''fct to plot station event data and all other stations within +-60min'''

    colrs_dict = pltcolor()  # create dictionary for colors

    # get all events, extracted from script _00_
    dfs_data_lst = list_all_full_path('.csv', path_to_events)

    shp_de = shapefile.Reader(shp_de)

    for i, event in enumerate(dfs_data_lst):

        (stn_one_id, ppt_stn_one, stns_2_ids,
         stn_one_xcoords, stn_one_ycoords,
         stns_2_xcoords, stns_2_ycoords,
         event_date, stns_2_ids_vals_dict,
         ppt_thr) = get_events_stn_data(event,
                                        path_stns_coords,
                                        xcoords_name,
                                        ycoords_name)

        coords_tuples = np.array([(x2, y2) for x2, y2
                                  in zip(stns_2_xcoords,
                                         stns_2_ycoords)])

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

        save_event_time = event_date.replace(
            ':', '_').replace(' ', '_').replace('-', '_')

        if not os.path.exists(os.path.join(path_to_dfs_simultaneous_events,
                                           'station_%s_ppt_thr_%smm_at_%s_2.png'
                                           % (str(int(stn_one_id)),
                                              ppt_thr, save_event_time))):
            print('plotting')
#
            fig, axes = plt.subplots(nrows=5, ncols=5,
                                     figsize=(80, 40), dpi=300,
                                     sharex=True, sharey=True)
#             axes = axes.ravel()

            ppt_stn_one_float_format = '% 0.2f' % ppt_stn_one
            texts = []
            for i in range(5):
                for j in range(5):
                    for shape_ in shp_de.shapeRecords():
                        lon = [i[0] for i in shape_.shape.points[:][::-1]]
                        lat = [i[1] for i in shape_.shape.points[:][::-1]]
#
                        x0, y0 = convert_coords_fr_wgs84_to_utm32_(
                            wgs82, utm32, lon, lat)
#
                        axes[i, j].scatter(x0,
                                           y0,
                                           marker='.',
                                           c='lightgrey',
                                           alpha=0.25,
                                           s=1)

                for i, stn2_id in enumerate(stns2_ids_in_circle):
                    stn2_xcoord = stns2_xcoords_in_circle[i]
                    stn2_ycoord = stns2_ycoords_in_circle[i]
                    time_idx_dict = stns_2_ids_vals_dict_in_circle[stn2_id]
#
                    for _, (idx, val) in enumerate(time_idx_dict.items()):
                        #
                        if len(val) > 0 and val[0] > 0:
                            #                                     # TODO fix Me
                            if idx == -60:
                                axes[0, 0].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

#                                 texts.append(axes[0, 0].text(stn2_xcoord,
#                                                              stn2_ycoord,
#                                                              val_float_format,
#                                                              color=colrs_dict[idx][0]))
                                adjust_text([axes[0, 0].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[0, 0])
                            if idx == -55:
                                axes[0, 1].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[0, 1].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[0, 1])
                            if idx == -50:
                                axes[0, 2].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[0, 2].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[0, 2])
                            if idx == -45:
                                axes[0, 3].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[0, 3].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[0, 3])
                            if idx == -40:
                                axes[0, 4].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])
#                                 print(val_float_format)
                                adjust_text([axes[0, 4].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[0, 4])
                            if idx == -35:
                                axes[1, 0].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))

                                val_float_format = '% 0.1f' % (val[0])
#                                 print(val_float_format, texts)
                                adjust_text([axes[1, 0].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[1, 0])

                            if idx == -30:
                                axes[1, 1].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[1, 1].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[1, 1])
                            if idx == -25:
                                axes[1, 2].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[1, 2].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[1, 2])

                            if idx == -20:
                                axes[1, 3].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[1, 3].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[1, 3])
                            if idx == -15:
                                axes[1, 4].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[1, 4].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[1, 4])
                            if idx == -10:
                                axes[2, 0].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[2, 0].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[2, 0])
                            if idx == -5:
                                axes[2, 1].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[2, 1].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[2, 1])
                            if idx == 0:
                                axes[2, 2].scatter(stn_one_xcoords,
                                                   stn_one_ycoords,
                                                   c='red',
                                                   marker='X',
                                                   s=25,
                                                   label=('Stn %s Ppt %0.2f mm'
                                                          % (str(int(stn_one_id)),
                                                             ppt_stn_one)))

                                axes[2, 2].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[2, 2].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0]),
                                    axes[2, 2].text(stn_one_xcoords,
                                                    stn_one_ycoords, ppt_stn_one_float_format,
                                                    color='k')], ax=axes[2, 2])
                            if idx == 5:
                                axes[2, 3].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[2, 3].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[2, 3])
                            if idx == 10:
                                axes[2, 4].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[2, 4].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[2, 4])
                            if idx == 15:
                                axes[3, 0].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[3, 0].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[3, 0])

                            if idx == 20:
                                axes[3, 1].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[3, 1].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[3, 1])
                            if idx == 25:
                                axes[3, 2].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[3, 2].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[3, 2])
                            if idx == 30:
                                axes[3, 3].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[3, 3].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[3, 3])
                            if idx == 35:
                                axes[3, 4].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[3, 4].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[3, 4])
                            if idx == 40:
                                axes[4, 0].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[4, 0].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[4, 0])
                            if idx == 45:
                                axes[4, 1].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[4, 1].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[4, 1])
                            if idx == 50:
                                axes[4, 2].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[4, 2].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[4, 2])
                            if idx == 55:
                                axes[4, 3].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[4, 3].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[4, 3])
                            if idx == 60:
                                axes[4, 4].scatter(stn2_xcoord, stn2_ycoord,
                                                   c=colrs_dict[idx][0],
                                                   marker=markers_time_dict[idx], s=15,
                                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                                          % (str(int(stn2_id)), val[0], idx)))
                                val_float_format = '% 0.1f' % (val[0])

                                adjust_text([axes[4, 4].text(
                                    stn2_xcoord,
                                    stn2_ycoord,
                                    val_float_format,
                                    color=colrs_dict[idx][0])], ax=axes[4, 4])
#                                     val_float_format = '% 0.1f mm %0.0f min' % (
#                                         val[0], idx)
#                                     val_float_format = '% 0.1f' % (val[0])
#
#                                     texts.append(ax.text(stn2_xcoord,
#                                                          stn2_ycoord,
#                                                          val_float_format,
#                                                          color=colrs_dict[idx][0]))
# #                         break
# #                 break

#
#             plt.set_title('Event_at_station_%s_ppt_thr_%smm_at_%s_30km_radius'
#                          % (str(int(stn_one_id)), ppt_thr, event_date))
#
#             ax.grid(alpha=0.25)
#
#             ax.set_xlabel("Longitude")
#             ax.set_ylabel("Latitude")
#
#             ax.set_xlim([250000, 950000])
#             ax.set_ylim([5200000, 6110000])
#             ax.set_aspect(1.0)

    #         ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),
    #                   fancybox=True, shadow=True, ncol=3)
#
#             plt.tight_layout()
            plt.savefig(os.path.join(path_to_dfs_simultaneous_events,
                                     'station_%s_ppt_thr_%smm_at_%s_2.png'
                                     % (str(int(stn_one_id)),
                                         ppt_thr, save_event_time)),
                        frameon=True, papertype='a4',
                        bbox_inches='tight', pad_inches=.2)
            plt.close()
#             break
        else:
            print('plots already exists')
            continue

        return


def plot_stations_in_convex_hull2(path_to_events,
                                  path_stns_coords,
                                  xcoords_name,
                                  ycoords_name,
                                  shp_de,
                                  radius=30000):
    '''fct to plot station event data and all other stations within +-60min'''

    colrs_dict = pltcolor()  # create dictionary for colors

    # get all events, extracted from script _00_
    dfs_data_lst = list_all_full_path('.csv', path_to_events)

    shp_de = shapefile.Reader(shp_de)

    for i, event in enumerate(dfs_data_lst):

        (stn_one_id, ppt_stn_one, stns_2_ids,
         stn_one_xcoords, stn_one_ycoords,
         stns_2_xcoords, stns_2_ycoords,
         event_date, stns_2_ids_vals_dict,
         ppt_thr) = get_events_stn_data(event,
                                        path_stns_coords,
                                        xcoords_name,
                                        ycoords_name)

        coords_tuples = np.array([(x2, y2) for x2, y2
                                  in zip(stns_2_xcoords,
                                         stns_2_ycoords)])

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

        save_event_time = event_date.replace(
            ':', '_').replace(' ', '_').replace('-', '_')

        time_vals = np.arange(-60, 61, 5)

        fig = plt.figure(figsize=(40, 20), dpi=75)
        # fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        # ax = fig.gca(projection='3d')
        ax = fig.add_subplot(111, projection='3d')
        ax.xaxis.pane.set_edgecolor('black')
        ax.yaxis.pane.set_edgecolor('black')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        ax.xaxis.set_major_locator(LinearLocator(10))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.04f'))
        ax.yaxis.set_major_locator(LinearLocator(10))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.04f'))

        for shape_ in shp_de.shapeRecords():
            lon = [i[0] for i in shape_.shape.points[:][::-1]]
            lat = [i[1] for i in shape_.shape.points[:][::-1]]
#
            x0, y0 = convert_coords_fr_wgs84_to_utm32_(
                wgs82, utm32, lon, lat)

            ax.scatter3D(x0, y0, 0, zdir='z',
                         color='k', alpha=0.05,
                         marker='.', s=25,  # linewidth=0.05,
                         label='')
#                         # plot first station

        ax.scatter3D(stn_one_xcoords,
                     stn_one_ycoords,
                     time_vals,
                     zdir='z',
                     c='k',
                     marker='.',
                     s=0)
        # label=('Stn %s Ppt %0.2f mm'
        #      % (str(int(stn_one_id)),
        #        ppt_stn_one)))
        ax.scatter3D(stn_one_xcoords,
                     stn_one_ycoords,
                     0,
                     zdir='z',
                     c='red',
                     marker='X',
                     s=25,
                     label=('Stn %s Ppt %0.2f mm'
                            % (str(int(stn_one_id)),
                               ppt_stn_one)))

        # plot all other simultaneous stations

        for i, stn2_id in enumerate(stns2_ids_in_circle):
            stn2_xcoord = stns2_xcoords_in_circle[i]
            stn2_ycoord = stns2_ycoords_in_circle[i]
            time_idx_dict = stns_2_ids_vals_dict_in_circle[stn2_id]
#
            for _, (idx, val) in enumerate(time_idx_dict.items()):
                #
                if len(val) > 0 and val[0] > 0:

                    x_vals = stn2_xcoord
                    y_vals = stn2_ycoord
                    z_vals = idx

                    ax.scatter3D(x_vals, y_vals, zs=z_vals, zdir='z',
                                 c=colrs_dict[idx][0], alpha=0.99,
                                 marker=markers_time_dict[idx], s=25,
                                 label=('Stn %s Ppt %0.2f mm'
                                        % (str(int(stn2_id)),
                                           val[0])))

        ax.zaxis.set_ticks(time_vals)
        ax.zaxis.set_ticklabels(time_vals)
        ax.set_xlabel('Longitude (x-axis)')
        ax.set_ylabel('Latitude (y-axis)')

        ax.set_title('Event_at_station_%s_ppt_thr_%smm_at_%s'
                     % (str(int(stn_one_id)), ppt_thr, event_date))

#         ax.set_xlim(min(lon), max(lon)), ax.set_ylim(min(lat), max(lat))
#         norm = mcolors.BoundaryNorm(ticks, cmap.N)
#         ax_legend = fig.add_axes([0.1725, 0.07525, 0.68, 0.0225], zorder=3)
#         cb = mpl.colorbar.ColorbarBase(ax_legend, ticks=ticks,
#                                        boundaries=ticks, norm=norm, cmap=cmap,
#                                        orientation='horizontal')
#         cb.set_label('Time_of_day_h', rotation=0)
        ax.set_aspect('auto')
#         cb.ax.set_xticklabels([str(i) for i in ticks])

        ax.view_init(25, 225)

#         ax.legend(loc=0)
#         cb.draw_all()
#         cb.set_alpha(1)
#         plt.tight_layout()
        plt.savefig(os.path.join(path_to_dfs_simultaneous_events,
                                 'station_%s_ppt_thr_%smm_at_%s_2_stns_in_%d.png'
                                 % (str(int(stn_one_id)),
                                     ppt_thr, save_event_time, radius)),
                    frameon=True, papertype='a4',
                    bbox_inches='tight', pad_inches=.2)
        plt.close()
        break
    return

# %%


def plot_stations_in_convex_hull3(path_to_events,
                                  path_stns_coords,
                                  xcoords_name,
                                  ycoords_name,
                                  shp_de,
                                  radius=30000):
    '''fct to plot station event data and all other stations within +-60min'''

    colrs_dict = pltcolor()  # create dictionary for colors

    # get all events, extracted from script _00_
    dfs_data_lst = list_all_full_path('.csv', path_to_events)

    shp_de = shapefile.Reader(shp_de)

    for i, event in enumerate(dfs_data_lst):

        (stn_one_id, ppt_stn_one, stns_2_ids,
         stn_one_xcoords, stn_one_ycoords,
         stns_2_xcoords, stns_2_ycoords,
         event_date, stns_2_ids_vals_dict,
         ppt_thr) = get_events_stn_data(event,
                                        path_stns_coords,
                                        xcoords_name,
                                        ycoords_name)

        coords_tuples = np.array([(x2, y2) for x2, y2
                                  in zip(stns_2_xcoords,
                                         stns_2_ycoords)])

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

        save_event_time = event_date.replace(
            ':', '_').replace(' ', '_').replace('-', '_')

        time_vals = np.arange(-60, 61, 5)

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

            for i, stn2_id in enumerate(stns2_ids_in_circle):
                stn2_xcoord = stns2_xcoords_in_circle[i]
                stn2_ycoord = stns2_ycoords_in_circle[i]
                time_idx_val = stns_2_ids_vals_dict_in_circle[stn2_id][timeval]

                if len(time_idx_val) > 0:

                    # plot first station
                    ax.scatter(stn_one_xcoords,
                               stn_one_ycoords, c='red',
                               marker='X',
                               s=50,
                               label=('Stn %s Ppt %0.2f mm'
                                      % (str(int(stn_one_id)),
                                         ppt_stn_one)))

                    ax.scatter(stn2_xcoord,
                               stn2_ycoord,
                               c=colrs_dict[timeval][0],
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
                         % (str(int(stn_one_id)), ppt_thr, event_date, timeval))

            ax.grid(alpha=0.25)

            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

            ax.set_xlim([250000, 950000])
            ax.set_ylim([5200000, 6110000])
            ax.set_aspect(1.0)

            plt.tight_layout()
            plt.savefig(os.path.join(path_to_dfs_simultaneous_events,
                                     'station_%s_ppt_thr_%smm_at_%s_2_%dmin_.png'
                                     % (str(int(stn_one_id)),
                                         ppt_thr, save_event_time,
                                         timeval)),
                        frameon=True, papertype='a4',
                        bbox_inches='tight', pad_inches=.2)
            plt.close()

    return


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    plot_stations_in_convex_hull3(path_to_dfs_simultaneous_events,
                                  path_to_ppt_coords_data, xcoord_name,
                                  ycoord_name, path_to_shpfile,
                                  30000)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
