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

from scipy.stats import pearsonr
plt.ioff()
#==============================================================================
#
#==============================================================================
main_dir = Path(os.getcwd())
os.chdir(main_dir)

# input dwd and netatmo interpolated rainfall fields

in_netatmo_ppt_nc_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\kriging_ppt_netatmo\netatmo_precipitation_kriging_2015-01-01_to_2019-06-01_1km_test_2.nc")

in_dwd_ppt_nc_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\kriging_ppt_netatmo\dwd_dwd_precipitation_kriging_2015-01-01_to_2019-06-01_1km_test_2.nc")


def read_nc_ppt_data(nc_file):
    ppt_dataset = nc.Dataset(nc_file)
    x_coords = ppt_dataset.variables['X'][:].data
    y_coords = ppt_dataset.variables['Y'][:].data

    ppt_data = ppt_dataset.variables['OK'][:].data

    time_data = ppt_dataset.variables['time'][:].data
    time_arr = nc.num2date(time_data, ppt_dataset.variables['time'].units,
                           calendar=ppt_dataset.variables['time'].calendar)
    time_ix = pd.DatetimeIndex(time_arr)
    return x_coords, y_coords, time_ix, ppt_data


def generate_correlation_map(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

    Returns
    -------
    np.array
      N X M array in which each element is a correlation coefficient.

    """
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must ' +
                         'have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    cov = np.dot(x,
                 y.T) - n * np.dot(mu_x[:, np.newaxis],
                                   mu_y[np.newaxis, :])
    return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    xnetatmo, ynetamo, time_netatmo, ppt_netatmo = read_nc_ppt_data(
        in_netatmo_ppt_nc_file)
    xdwd, ydwd, time_dwd, ppt_dwd = read_nc_ppt_data(in_dwd_ppt_nc_file)

    assert all(xnetatmo) == all(xdwd)
    assert all(ynetamo) == all(ydwd)
    diff_netatmo_dwd_ppt = ppt_netatmo - ppt_dwd
    divide_netatmo_dwd_ppt = np.divide(ppt_netatmo, ppt_dwd)

    # for corrrelation
    xarr = np.linspace(1, len(ydwd), len(ydwd))
    yarr = np.linspace(1, len(ydwd), len(ydwd))

    for ix in range(len(time_dwd)):
        print(time_dwd[ix])
        fig = plt.figure(figsize=(12, 8), dpi=150)
        plt.pcolormesh(xnetatmo, ynetamo,
                       diff_netatmo_dwd_ppt[ix],
                       cmap=plt.get_cmap('jet'),
                       vmin=-100,
                       vmax=100)
        plt.title('Netatmo - DWD (ppt/m) for %s' % str(time_dwd[ix]))
        plt.colorbar(extend='both', label='ppt/m')
        plt.grid(alpha=0.25)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\compare_kriging_difference_and_correlation_maps\netatmo_diff_dwd_%s_.png"
                    % str(time_dwd[ix]).replace(':', '_').replace(' ', '_'))
        plt.close('all')
        #======================================================================
        #
        #======================================================================
        fig = plt.figure(figsize=(12, 8), dpi=150)
        plt.pcolormesh(xnetatmo, ynetamo,
                       divide_netatmo_dwd_ppt[ix],
                       cmap=plt.get_cmap('jet'),
                       vmin=-100,
                       vmax=100)
        plt.title('Netatmo / DWD (ppt/m) for %s' % str(time_dwd[ix]))
        plt.colorbar(extend='both', label='ppt/m')
        plt.grid(alpha=0.25)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\compare_kriging_difference_and_correlation_maps\netatmo_div_dwd_%s_.png"
                    % str(time_dwd[ix]).replace(':', '_').replace(' ', '_'))
        plt.close('all')

        #======================================================================
        #
        #======================================================================

        fig = plt.figure(figsize=(12, 8), dpi=150)
        plt.pcolormesh(xarr, yarr,
                       generate_correlation_map(ppt_netatmo[ix], ppt_dwd[ix]),
                       cmap=plt.get_cmap('jet'),
                       vmin=-1,
                       vmax=1)
        plt.title('Correlation Netatmo DWD (ppt/m) for %s' % str(time_dwd[ix]))
        plt.colorbar(extend='both', label='ppt/m')
        plt.grid(alpha=0.25)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\compare_kriging_difference_and_correlation_maps\netatmo_corr_dwd_%s_.png"
                    % str(time_dwd[ix]).replace(':', '_').replace(' ', '_'))
        plt.close('all')

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
