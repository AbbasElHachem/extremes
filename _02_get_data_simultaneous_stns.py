# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
GOAL: get station ids and coordinates for plotting on map

Get for every event the 'main' station and all other 
station where rainfall was above the threshold, occurring
with in +-60min, return the result to be used in script _03_
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# ===================================================

from _01_compare_pairwise_stations import time_shifts_arr_floats

import os

import timeit
import time

import fnmatch
import pandas as pd

#==============================================================================
#
#==============================================================================


def list_all_full_path(ext, file_dir):
    """
    Purpose: To return full path of files in all dirs of a 
            given folder with a
            given extension in ascending order.
    Description of the arguments:
        ext (string) = Extension of the files to list
            e.g. '.txt', '.tif'.
        file_dir (string) = Full path of the folder in which the files
            reside.
    """
    new_list = []
    patt = '*' + ext
    for root, _, files in os.walk(file_dir):
        for elm in files:
            if fnmatch.fnmatch(elm, patt):
                full_path = os.path.join(root, elm)
                new_list.append(full_path)
    return(sorted(new_list))


#==============================================================================
#
#==============================================================================


def read_df_stns_coord(path_to_df):
    ''' read df_file for stations ids and coordinates'''
    in_df_coords = pd.read_csv(path_to_df, sep=';',
                               index_col=3, engine='c')
    stns_ids = in_df_coords.index.values

    return stns_ids, in_df_coords

#==============================================================================
#
#==============================================================================


def get_events_stn_data(df_event_file,
                        path_stns_coords,
                        xcoord_name,
                        ycoord_name):
    ''' fct to get simultaneous evetns, at station one and all other stns'''
    stn_ids, df_coords = read_df_stns_coord(path_stns_coords)
    print(df_event_file)
    split_file = df_event_file.split('_')

    stn_one_id = float(split_file[-3])
    assert stn_one_id in stn_ids

    ppt_thr = split_file[-1].replace('.csv', '')
    # print(split_file)
    ppt_stn_one = float(split_file[3])

    stn_one_xcoords = df_coords.loc[stn_one_id, xcoord_name]
    stn_one_ycoords = df_coords.loc[stn_one_id, ycoord_name]

    event_date = (split_file[-8] + ' ' + split_file[-7] + ':' +
                  split_file[-6] + ':' + split_file[-5])

    in_df_event = pd.read_csv(df_event_file,
                              sep=',',
                              index_col=0,
                              engine='c')
    stns_2_ids = [float(col) for col in in_df_event.columns.values]
    print(in_df_event)
    # assert stns_2_ids in stn_ids

    stns_2_ids_vals_dict = {stn2:
                            {time_shift: []
                             for time_shift in time_shifts_arr_floats}
                            for stn2 in stns_2_ids}

    stns_2_xcoords, stns_2_ycoords = [], []

    for stn2_id in stns_2_ids:

        stns_2_xcoords.append(df_coords.loc[stn2_id, xcoord_name])
        stns_2_ycoords.append(df_coords.loc[stn2_id, ycoord_name])

        time_idx_shift = in_df_event.loc[:,
                                         str(int(stn2_id))].dropna().index
        for shift in time_idx_shift:
            stns_2_ids_vals_dict[stn2_id][shift].append(
                in_df_event.loc[shift, str(int(stn2_id))])

    return (stn_one_id, ppt_stn_one, stns_2_ids,
            stn_one_xcoords, stn_one_ycoords,
            stns_2_xcoords, stns_2_ycoords,
            event_date, stns_2_ids_vals_dict, ppt_thr)


#==============================================================================
#
#==============================================================================


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
