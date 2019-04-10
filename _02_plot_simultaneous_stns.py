# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
GOAL: get station ids and coordinates for plotting on map

"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# ===================================================


import timeit
import time

import pandas as pd

path_to_ppt_coords_data = (r'X:\exchange\ElHachem'
                           r'\niederschlag_deutschland'
                           r'\1993_2016_5min_merge_nan.csv')


def read_df_stns_coord(path_to_df):
    in_df_coords = pd.read_csv(path_to_df, sep=';')
    stns_ids = in_df_coords['Stationsnr.'].values
    stns_xcoords = in_df_coords['Rechtswert'].values
    stns_ycoords = in_df_coords['Hochwert'].values
    return stns_ids, stns_xcoords, stns_ycoords


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program
    stns_ids, stns_xcoords, stns_ycoords = read_df_stns_coord(
        path_to_ppt_coords_data)
    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
