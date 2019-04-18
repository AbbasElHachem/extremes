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


import shapefile
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as spatial

from adjustText import adjust_text

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
            fig = plt.figure(figsize=(15, 15), dpi=150)

            ax = fig.add_subplot(111)

            for shape_ in shp_de.shapeRecords():
                lon = [i[0] for i in shape_.shape.points[:][::-1]]
                lat = [i[1] for i in shape_.shape.points[:][::-1]]

                x0, y0 = convert_coords_fr_wgs84_to_utm32_(
                    wgs82, utm32, lon, lat)

                ax.scatter(x0, y0, marker='.', c='lightgrey', alpha=0.5, s=2)
            # plot first station
            ax.scatter(stn_one_xcoords, stn_one_ycoords, c='red',
                       marker='X', s=50,
                       label=('Stn %s Ppt %0.2f mm'
                              % (str(int(stn_one_id)),
                                 ppt_stn_one)))

            # plot all other simultaneous stations
            texts = []
            for i, stn2_id in enumerate(stns2_ids_in_circle):
                stn2_xcoord = stns2_xcoords_in_circle[i]
                stn2_ycoord = stns2_ycoords_in_circle[i]
                time_idx_dict = stns_2_ids_vals_dict_in_circle[stn2_id]

                for _, (idx, val) in enumerate(time_idx_dict.items()):

                    if len(val) > 0 and val[0] > 0:
                        # TODO fix Me

                        ax.scatter(stn2_xcoord, stn2_ycoord, c=colrs_dict[idx][0],
                                   marker=markers_time_dict[idx], s=45,
                                   label=('Stn %s Ppt %0.2f at %0.0f min'
                                          % (str(int(stn2_id)), val[0], idx)))

                        val_float_format = '% 0.1f mL %0.0f min' % (
                            val[0], idx)
    #                     val_float_format = '% 0.1f' % (val[0])

                        texts.append(ax.text(stn2_xcoord,
                                             stn2_ycoord,
                                             val_float_format,
                                             color=colrs_dict[idx][0]))

            ppt_stn_one_float_format = '% 0.2f' % ppt_stn_one
            texts.append(ax.text(stn_one_xcoords,
                                 stn_one_ycoords,
                                 ppt_stn_one_float_format))

            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle='->', color='red', lw=0.25))

            ax.set_title('Event_at_station_%s_ppt_thr_%smm_at_%s_30km_radius'
                         % (str(int(stn_one_id)), ppt_thr, event_date))

            ax.grid(alpha=0.25)

            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

            ax.set_xlim([250000, 950000])
            ax.set_ylim([5200000, 6110000])
            ax.set_aspect(1.0)

    #         ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),
    #                   fancybox=True, shadow=True, ncol=3)

            plt.tight_layout()
            plt.savefig(os.path.join(path_to_dfs_simultaneous_events,
                                     'station_%s_ppt_thr_%smm_at_%s_2.png'
                                     % (str(int(stn_one_id)),
                                         ppt_thr, save_event_time)),
                        frameon=True, papertype='a4',
                        bbox_inches='tight', pad_inches=.2)
            plt.close()
            break
        else:
            print('plots already exists')
            continue

        return


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    plot_stations_in_convex_hull(path_to_dfs_simultaneous_events,
                                 path_to_ppt_coords_data, xcoord_name,
                                 ycoord_name, path_to_shpfile,
                                 30000)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
