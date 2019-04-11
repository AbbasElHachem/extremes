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

# ===================================================

from _02_plot_simultaneous_stns import list_all_full_path
from _02_plot_simultaneous_stns import get_events_stn_data


import os
import timeit
import time

import shapefile
import pyproj
import matplotlib.pyplot as plt

import numpy as np

from matplotlib import colors as mcolors

#==============================================================================
#
#==============================================================================
path_to_dfs_simultaneous_events = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'

path_to_ppt_coords_data = (r'X:\exchange\ElHachem'
                           r'\niederschlag_deutschland'
                           r'\1993_2016_5min_merge_nan.csv')

path_to_shpfile = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                   r'\shp_file_germany\DEU_adm1.shp')
xcoord_name = 'Rechtswert'
ycoord_name = 'Hochwert'

# def epsg wgs84 and utm32
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

#==============================================================================
# Get colors by name
#==============================================================================

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Sort colors by hue, saturation, value and name.
by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                for name, color in colors.items())
sorted_names_all = [name for hsv, name in by_hsv]


markers = ['o', '.', ',', 'x', '+', 'v', '^', '<', '>', 's', 'd',
           '1', '8', '8', '1', 'd', 's', '>', '<', '^', 'v', '+',
           'x', ',', '.', 'o']

markers_time_dict = {i: [] for i in np.arange(-60, 61, 5)}
for i, m in zip(np.arange(-60, 61, 5), markers):
    markers_time_dict[i] = m

#==============================================================================
#
#==============================================================================


def convert_coords_fr_wgs84_to_utm32_(epgs_initial_str, epsg_final_str,
                                      first_coord, second_coord):
    '''fct to convert points from wgs 84 to utm32'''
    initial_epsg = pyproj.Proj(epgs_initial_str)
    final_epsg = pyproj.Proj(epsg_final_str)
    x, y = pyproj.transform(initial_epsg, final_epsg,
                            first_coord, second_coord)
    return x, y

#==============================================================================
#
#==============================================================================

# Function to map the colors as a list from the input list of x variables


def pltcolor():
    cols_time_dict = {i: [] for i in np.arange(-60, 61, 5)}

    for i, l in enumerate(cols_time_dict.keys()):
        if l == -60:
            cols_time_dict[l].append(sorted_names_all[i * 12])
        elif l == -55:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -50:
            cols_time_dict[l].append(sorted_names_all[i * 8])
        elif l == -45:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -40:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -35:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -30:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -25:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -20:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -15:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -10:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -5:
            cols_time_dict[l].append(sorted_names_all[i * 6])

        elif l == 0:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 5:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 10:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 15:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 20:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 25:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 30:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 35:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 40:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 45:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 50:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 55:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 60:
            cols_time_dict[l].append(sorted_names_all[i * 6])
    return cols_time_dict

#==============================================================================
#
#==============================================================================


def plot_coordinates(path_to_events,
                     path_stns_coords,
                     xcoords_name,
                     ycoords_name,
                     shp_de):

    colrs_dict = pltcolor()

    dfs_data_lst = list_all_full_path('.csv', path_to_events)

    shp_de = shapefile.Reader(shp_de)

    for i, event in enumerate(dfs_data_lst):

        (stn_one_id, ppt_stn_one, stns_2_ids,
         stn_one_xcoords, stn_one_ycoords,
         stns_2_xcoords, stns_2_ycoords,
         event_date, stns_2_ids_vals_dict) = get_events_stn_data(event,
                                                                 path_stns_coords,
                                                                 xcoords_name,
                                                                 ycoords_name)

        fig = plt.figure(figsize=(20, 12), dpi=150)

        ax = fig.add_subplot(111)
        for shape_ in shp_de.shapeRecords():
            lon = [i[0] for i in shape_.shape.points[:][::-1]]
            lat = [i[1] for i in shape_.shape.points[:][::-1]]

            x0, y0 = convert_coords_fr_wgs84_to_utm32_(wgs82, utm32, lon, lat)

            ax.scatter(x0, y0, marker='.', c='grey', alpha=0.35, s=2)

        ax.scatter(stn_one_xcoords, stn_one_ycoords, c='r', marker='x', s=40,
                   label=('Stn %s ppt was %0.2f mm'
                          % (str(int(stn_one_id)), ppt_stn_one)))

        for i, stn2_id in enumerate(stns_2_ids):
            stn2_xcoord = stns_2_xcoords[i]
            stn2_ycoord = stns_2_ycoords[i]
            time_idx_dict = stns_2_ids_vals_dict[stn2_id]

            for _, (idx, val) in enumerate(time_idx_dict.items()):

                if len(val) > 0:
                    ax.scatter(stn2_xcoord, stn2_ycoord, c=colrs_dict[idx],
                               marker=markers_time_dict[idx], s=40,
                               label=('Stn %s ppt was %0.2f at %0.0f min'
                                      % (str(int(stn2_id)), val[0], idx)))

        ax.set_title(event_date)
        ax.legend(loc=0)
        plt.tight_layout()

        plt.savefig(os.path.join(path_to_dfs_simultaneous_events,
                                 '_event_%d_.png' % i),
                    frameon=True, papertype='a4',
                    bbox_inches='tight', pad_inches=.2)
        plt.close()
#         break
#         return


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program
    colrs_dict = pltcolor()
    plot_coordinates(path_to_dfs_simultaneous_events,
                     path_to_ppt_coords_data, xcoord_name,
                     ycoord_name, path_to_shpfile)
    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
