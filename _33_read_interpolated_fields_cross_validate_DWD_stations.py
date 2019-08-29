# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Compare interpolated rainfall fields
Purpose: Calculate Difference and Correlation Maps between two interpolations

Created on: 2019-08-08

Parameters
----------
Interpolated DWD rainfall data
Interpolated Netatmo rainfall data

Returns
-------
Difference maps (minus, divide) for every cell
Correlation map
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

from pathlib import Path

import os
import timeit
import time

import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd
import numpy as np

from _31_read_interpolated_ppt_fields_and_calc_diff_map_corr_map import read_nc_ppt_data
plt.ioff()
#==============================================================================
#
#==============================================================================
main_dir = Path(os.getcwd())
os.chdir(main_dir)

# input dwd and netatmo interpolated rainfall fields

in_netatmo_netatmo_ppt_nc_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\kriging_ppt_netatmo"
    r"\Netatmo_dwd_daily_precipitation_kriging_2015-01-01_to_2019-06-01_1km_mid_rg_gd_stns.nc")

# DWD stations coordinates UTM32

in_stn_coords_df_loc = os.path.join(
    r"F:\download_DWD_data_recent"
    r"\station_coordinates_names_hourly_only_in_BW_utm32.csv")

if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    xnetatmo, ynetamo, time_netatmo, ppt_netatmo = read_nc_ppt_data(
        in_netatmo_netatmo_ppt_nc_file)

    in_stns_coords_df = pd.read_csv(
        in_stn_coords_df_loc,
        sep=';',
        index_col=0,
        encoding='utf-8')

    xdwd = in_stns_coords_df.loc[:, 'X'].values.ravel()
    ydwd = in_stns_coords_df.loc[:, 'Y'].values.ravel()

    xdiff_lst = xnetatmo - xdwd[0]
#     np.where(np.isclose(xdwd[0], xnetatmo))
    plt.ioff()
    plt.plot(xnetatmo)
    plt.plot(xdwd)
    plt.scatter(xdwd, ydwd, c='r', marker='+', alpha=0.75)
    plt.scatter(xnetatmo, ynetamo, c='b', marker=',', alpha=0.25)
    plt.show()
    print('Cross validating DWD data')
# #     assert all(xnetatmo) == all(xdwd)
#     assert all(ynetamo) == all(ydwd)
# #     diff_netatmo_dwd_ppt = ppt_netatmo - ppt_dwd
#     divide_netatmo_dwd_ppt = np.divide(ppt_netatmo, ppt_dwd)
#
#     # for corrrelation
#     xarr = np.linspace(1, len(ydwd), len(ydwd))
#     yarr = np.linspace(1, len(ydwd), len(ydwd))
#
#     for ix in range(len(time_dwd)):
#
#         print(time_dwd[ix])
#         #======================================================================
#         # difference
#         #======================================================================
#         fig = plt.figure(figsize=(12, 8), dpi=150)
#         plt.pcolormesh(xnetatmo, ynetamo,
#                        diff_netatmo_dwd_ppt[ix],
#                        cmap=plt.get_cmap('jet'),
#                        vmin=-100,
#                        vmax=100)
#         plt.title('Netatmo - DWD (ppt/m) for %s' % str(time_dwd[ix]))
#         plt.colorbar(extend='both', label='ppt/m')
#         plt.grid(alpha=0.25)
#         plt.xlabel('X')
#         plt.ylabel('Y')
#         plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\compare_kriging_difference_and_correlation_maps\netatmo_diff_dwd_%s_.png"
#                     % str(time_dwd[ix]).replace(':', '_').replace(' ', '_'))
#         plt.close('all')
#         #======================================================================
#         # Divide
#         #======================================================================
#         fig = plt.figure(figsize=(12, 8), dpi=150)
#         plt.pcolormesh(xnetatmo, ynetamo,
#                        divide_netatmo_dwd_ppt[ix],
#                        cmap=plt.get_cmap('jet'),
#                        vmin=-10,
#                        vmax=10)
#         plt.title('Netatmo / DWD (ppt/m) for %s' % str(time_dwd[ix]))
#         plt.colorbar(extend='both', label='ppt/m')
#         plt.grid(alpha=0.25)
#         plt.xlabel('X')
#         plt.ylabel('Y')
#         plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\compare_kriging_difference_and_correlation_maps\netatmo_div_dwd_%s_.png"
#                     % str(time_dwd[ix]).replace(':', '_').replace(' ', '_'))
#         plt.close('all')

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
