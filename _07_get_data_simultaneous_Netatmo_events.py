# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
GOAL: get Netatmo station ids and coordinates for plotting on map

Get for every event the 'main' station and all other 
station where rainfall was above the threshold, occurring
with in +-120min, return the result to be used in script _05_
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

import pandas as pd

from _00_additional_functions import get_stns_ids_coords_file
from _06_NetatmoData_get_simultaneous_events import time_shifts_arr_floats


#==============================================================================
#
#==============================================================================


def read_df_stns_coord(path_to_df):
    ''' read df_file for stations ids and coordinates'''
    in_df_coords = pd.read_csv(path_to_df, sep=';',
                               index_col=0, engine='c')
    return in_df_coords


#==============================================================================
#
#==============================================================================


def get_netatmo_events_stn_data(df_event_file,
                                path_stns_coords,
                                xcoord_name,
                                ycoord_name):
    ''' fct to get simultaneous evetns, at station one and all other stns'''
    df_coords = read_df_stns_coord(path_stns_coords)
    stn_ids = get_stns_ids_coords_file(path_stns_coords)

    assert len(stn_ids) > 0, 'error in gettiog all stn ids'

    print(df_event_file)
    split_file = df_event_file.split('_')

    stn_one_id = df_event_file[-28:-11]
    print('First station ID is', stn_one_id)
    assert stn_one_id in stn_ids, 'wrong stn id'

    ppt_thr = split_file[-1].replace('.csv', '')
    # print(split_file)
    ppt_stn_one = float(split_file[3])
    stn_one_id_orig = stn_one_id.replace('_', ':')

    assert xcoord_name in df_coords.columns, 'wrong column for coords'

    stn_one_xcoords = df_coords.loc[stn_one_id_orig, xcoord_name]
    stn_one_ycoords = df_coords.loc[stn_one_id_orig, ycoord_name]

    event_date = (split_file[6] + ' ' + split_file[7] + ':' +
                  split_file[8] + ':' + split_file[9])

    in_df_event = pd.read_csv(df_event_file,
                              sep=',',
                              index_col=0,
                              engine='c')

    stns_2_ids = [str(col) for col in in_df_event.columns.values
                  if col.replace('_', ':') in df_coords.index]
    print(stns_2_ids)
    # assert stns_2_ids in stn_ids

    stns_2_ids_vals_dict = {stn2:
                            {time_shift: []
                             for time_shift in time_shifts_arr_floats}
                            for stn2 in stns_2_ids}

    stns_2_xcoords, stns_2_ycoords = [], []

    for stn2_id in stns_2_ids:
        stn2_id_orig = stn2_id.replace('_', ':')
        if stn2_id_orig in df_coords.index:
            stns_2_xcoords.append(df_coords.loc[stn2_id_orig, xcoord_name])
            stns_2_ycoords.append(df_coords.loc[stn2_id_orig, ycoord_name])

            time_idx_shift = in_df_event.loc[:,
                                             str(stn2_id)].dropna().index
            for shift in time_idx_shift:
                stns_2_ids_vals_dict[stn2_id][shift].append(
                    in_df_event.loc[shift, str(stn2_id)])

    return (stn_one_id, ppt_stn_one, stns_2_ids,
            stn_one_xcoords, stn_one_ycoords,
            stns_2_xcoords, stns_2_ycoords,
            event_date, stns_2_ids_vals_dict, ppt_thr)
#==============================================================================
#
#==============================================================================
