# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Inverse Distance Weighting
Purpose: Module for the Inverse Distance Weighting interpolation technique

Created on: 2019-07-22

The main idea of this interpolation strategy lies in fact that it is not
desirable to honour local high/low values but rather to look at a moving
average of nearby data points and estimate the local trends. 
The node value is calculated by averaging the weighted sum of all the points.
Data points that lie progressively farther from the node inuence much less
the computed value than those lying closer to the node.

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

# =============================================================================

from pathlib import Path

import os
import timeit
import time
import pyproj
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _00_additional_functions import convert_coords_fr_wgs84_to_utm32_

main_dir = Path(os.getcwd())
os.chdir(main_dir)


#==============================================================================
#
#==============================================================================
path_to_ppt_netatmo_data = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_ppt_netatmo_data), 'wrong NETATMO Ppt file'

path_to_netatmo_coords_df_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\rain_bw_1hour"
    r"\netatmo_bw_1hour_coords.csv")
assert os.path.exists(path_to_netatmo_coords_df_file), 'wrong DWD coords file'


# used in Netatmo coords df
lon_col_name = ' lon'
lat_col_name = ' lat'

# def epsg wgs84 and utm32 for coordinates conversion
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

print('\n######\n reading all dfs \n#######\n')
# read netatmo ppt df
in_netatmo_ppt_stns_df = pd.read_csv(path_to_ppt_netatmo_data,
                                     index_col=0, sep=';',
                                     engine='c')
in_netatmo_ppt_stns_df.index = pd.to_datetime(in_netatmo_ppt_stns_df.index,
                                              format='%Y-%m-%d %H:%M:%S')
stns_ppt_ids = in_netatmo_ppt_stns_df.columns

stn_inter = stns_ppt_ids[10]

in_netatmo_ppt_stns_df_obs = in_netatmo_ppt_stns_df.loc[
    :, in_netatmo_ppt_stns_df.columns != stn_inter]


# read netatmo ppt coords df (for plotting)
in_netatmo_df_coords = pd.read_csv(path_to_netatmo_coords_df_file, sep=';',
                                   index_col=0, engine='c')

in_netatmo_df_coords.index = [stnid.replace(
    ':', '_') for stnid in in_netatmo_df_coords.index]

# TODO: Intersect two indexes
lon_stn_netatmo = in_netatmo_df_coords.loc[
    in_netatmo_df_coords.index == in_netatmo_ppt_stns_df_obs.columns,
    lon_col_name].values

lat_stn_netatmo = in_netatmo_df_coords.loc[
    in_netatmo_df_coords.index != stn_inter, lat_col_name].values

ppt_data = in_netatmo_ppt_stns_df_obs.loc['2018-08-30 00:00:00', :]

xstns, ystns = convert_coords_fr_wgs84_to_utm32_(wgs82, utm32,
                                                 lon_stn_netatmo,
                                                 lat_stn_netatmo)

lon_inter = in_netatmo_df_coords.loc[stn_inter, lon_col_name]
lat_inter = in_netatmo_df_coords.loc[stn_inter, lat_col_name]

x_inter, y_inter = convert_coords_fr_wgs84_to_utm32_(wgs82, utm32,
                                                     lon_inter,
                                                     lat_inter)
plt.ioff()
plt.scatter(x_inter, y_inter, c='r', marker='+')
plt.scatter(xstns, ystns, c='b', marker='.', alpha=.25)
plt.show()
#==============================================================================
#
#==============================================================================


def distance_matrix(x0, y0, x1, y1):
    obs = np.vstack((x0, y0)).T
    interp = np.vstack((x1, y1)).T

    # Make a distance matrix between pairwise observations
    # Note: from <http://stackoverflow.com/questions/1871536>
    # (Yay for ufuncs!)
    d0 = np.subtract.outer(obs[:, 0], interp[:, 0])
    d1 = np.subtract.outer(obs[:, 1], interp[:, 1])
    return np.hypot(d0, d1)


distances = distance_matrix(xstns, ystns,
                            xstns[10], ystns[10])


def simple_idw(x, y, z, xi, yi):
    dist = distance_matrix(x, y, xi, yi)

    # In IDW, weights are 1 / distance
    weights = 1.0 / dist

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Multiply the weights for each interpolated point
    # by all observed Z-values
    zi = np.dot(weights.T, z)
    return zi


simple_idw(xstns, ystns, ppt_data, x_inter, y_inter)


def plot(x, y, z, grid):
    plt.ioff()
    plt.figure()
    plt.imshow(grid, extent=(x.min(), x.max(), y.max(), y.min()))
    plt.scatter(x, y, c=z)
    plt.colorbar()


def main():
    # Setup: Generate data...
    n = 10
    nx, ny = 50, 50
    x, y, z = map(np.random.random, [n, n, n])
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    xi, yi = np.meshgrid(xi, yi)
    xi, yi = xi.flatten(), yi.flatten()

    # Calculate IDW
    grid1 = simple_idw(x, y, z, xi, yi)
    grid1 = grid1.reshape((ny, nx))
    plot(x, y, z, grid1)
    plt.title('IDW')
    plt.show()

#==============================================================================
#
#==============================================================================


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

#     main()

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
