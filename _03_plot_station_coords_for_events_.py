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

from itertools import cycle

from _02_plot_simultaneous_stns import list_all_full_path
from _02_plot_simultaneous_stns import get_events_stn_data


import os
import timeit
import time

import shapefile
import pyproj
import matplotlib.pyplot as plt
#import numpy as np

cycol = cycle(['blue', 'green', 'orange', 'cyan', 'm',
               'orange', 'darkblue', 'lime', 'gold', 'maroon',
               'olive', 'violet', 'tomato'])

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

def plot_coordinates(path_to_events,
                     path_stns_coords,
                     xcoords_name,
                     ycoords_name,
                     shp_de):
    dfs_data_lst = list_all_full_path('.csv', path_to_events)

    shp_de = shapefile.Reader(shp_de)

    for i, event in enumerate(dfs_data_lst):

        (stn_one_id, ppt_stn_one,
         stns_2_ids,
         stn_one_xcoords, stn_one_ycoords,
         stns_2_xcoords, stns_2_ycoords,
         event_date, stns_2_ids_vals_dict) = get_events_stn_data(event,
                                                                 path_stns_coords,
                                                                 xcoords_name,
                                                                 ycoords_name)

        fig = plt.figure(figsize=(20, 12), dpi=75)

        ax = fig.add_subplot(111)
        for shape_ in shp_de.shapeRecords():
            lon = [i[0] for i in shape_.shape.points[:][::-1]]
            lat = [i[1] for i in shape_.shape.points[:][::-1]]

            x0, y0 = convert_coords_fr_wgs84_to_utm32_(wgs82, utm32, lon, lat)

            ax.scatter(x0, y0, marker='.', c='grey', alpha=0.35, s=2)

        ax.scatter(stn_one_xcoords, stn_one_ycoords, c='r', marker='x', s=40,
                   label=('Stn %s ppt was %0.2f'
                          % (str(int(stn_one_id)), ppt_stn_one)))

        for i, stn2_id in enumerate(stns_2_ids):
            stn2_xcoord = stns_2_xcoords[i]
            stn2_ycoord = stns_2_ycoords[i]
            time_idx_dict = stns_2_ids_vals_dict[stn2_id]

            for _, (idx, val) in enumerate(time_idx_dict.items()):
                if len(val) > 0:
                    ax.scatter(stn2_xcoord, stn2_ycoord, c=next(cycol),
                               marker='8', s=40,
                               label=('Stn %s ppt was %0.2f at %0.0f min'
                                      % (str(int(stn2_id)), val[0], idx)))

        ax.set_title(event_date)
        ax.legend(loc=0)
        plt.savefig(os.path.join(path_to_dfs_simultaneous_events,
                                 '_event_%d_.png' % i))
#         return


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    plot_coordinates(path_to_dfs_simultaneous_events,
                     path_to_ppt_coords_data, xcoord_name, ycoord_name, path_to_shpfile)
    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
