# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
1) for every Netatmo precipitation station
select the closest Netatmo temperature station, 
intersect the two dataframes, select the days where only temperature
was above 1°c

2)from the new Netatmo ppt stations, readthe closeset DWD station
calcualte for the staions the values of P1 and compare the two
if above, below or equal (+-10 %),

3) plot the results on a map,
 see if the NEtatmo stations behave similar or not

"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

from pathlib import Path

import os
import timeit
import time

import numpy as np
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import matplotlib.dates as md

import shapefile

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears


from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from pandas.plotting import register_matplotlib_converters


from _00_additional_functions import (calculate_probab_ppt_below_thr,
                                      resample_intersect_2_dfs)

from b_get_data import HDF5

#==============================================================================
#
#==============================================================================
main_dir = Path(os.getcwd())
os.chdir(main_dir)

register_matplotlib_converters()

plt.ioff()

rc('font', size=13)
rc('font', family='serif')
rc('axes', labelsize=13)
rcParams['axes.labelpad'] = 13

majorLocator = MultipleLocator(2)
minorLocator = MultipleLocator(1)
#==============================================================================
#
#==============================================================================
path_to_ppt_netatmo_data = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_ppt_netatmo_data), 'wrong NETATMO Ppt file'

path_to_temp_netatmo_data = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\temperature_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_temp_netatmo_data), 'wrong NETATMO Temp file'


distance_matrix_df_file_ppt_temp = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\NetAtmo_BW\distance_mtx_in_m_NetAtmo_ppt_Netatmo_temp.csv")
assert os.path.exists(distance_matrix_df_file_ppt_temp), 'wrong distance df'


path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')
assert os.path.exists(path_to_ppt_hdf_data), 'wrong DWD Ppt file'

distance_matrix_netatmo_dwd_df_file = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD.csv')
assert os.path.exists(
    distance_matrix_netatmo_dwd_df_file), 'wrong Distance MTX  file'

path_to_netatmo_coords_df_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\rain_bw_1hour"
    r"\netatmo_bw_1hour_coords.csv")
assert os.path.exists(path_to_netatmo_coords_df_file), 'wrong DWD coords file'

dwd_coords_df_file = (r'X:\exchange\ElHachem\niederschlag_deutschland'
                      r'\1993_2016_5min_merge_nan.csv')
assert os.path.exists(dwd_coords_df_file), 'wrong DWD coords file'

path_to_shpfile = (
    r"X:\exchange\ElHachem\Netatmo\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp")

assert os.path.exists(path_to_shpfile), 'wrong shapefile path'

out_save_dir_orig = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

#==============================================================================
#
#==============================================================================
x_col_name = ' lon'
y_col_name = ' lat'


# threshold for max ppt value per hour
max_ppt_thr = 100.
ppt_min_thr = 1  # used when calculating p0

ppt_thrs_list = [0.5, 1, 2, 0.5, 5]

# till 1 day '5min', '10min', '15min', '30min',
aggregation_frequencies = ['60min', '90min', '120min', '180min', '240min',
                           '360min', '480min', '720min', '1440min']
# aggregation_frequencies = ['60min']

use_temp_thr = True
#==============================================================================
#
#==============================================================================


def compare_p1_dwd_p1_netatmo(val_dwd, val_netatmo):
    '''
    use this function to compare two values from two different stations
    find if the values are equal +-10% , if smaller or if bigger
    return the result for being saved in a dataframe
    plot the result as scatter plots

    marker, color 
    '''
    _10_percent_dwd = 0.1 * val_dwd

    if val_dwd - _10_percent_dwd <= val_netatmo <= val_dwd + _10_percent_dwd:
        return 's', 'b'  # '='
    if val_netatmo < val_dwd - _10_percent_dwd:
        return '_', 'g'  # '-'
    if val_dwd + _10_percent_dwd < val_netatmo:
        return '+', 'r'  # '+'
    return

#==============================================================================
#
#==============================================================================


def select_netatmo_ppt_abv_netatmo_temp_thr(
        netatmo_ppt_df_file,  # path to df of all netatmo ppt stations
        netatmo_temperature_df_file,  # path to df of all netatmo temp stations
        path_to_dwd_data,  # path to dwd ppt hdf5 data
        netatmo_ppt_coords_df,  # path to netatmo ppt coords df
        distance_matrix_netatmo_ppt_netatmo_temp,
        distance_matrix_netatmo_ppt_dwd_ppt,
        temp_freq_resample,
        use_temp_thr=use_temp_thr,  # remove all data where temp<temp_thr
        temp_thr=1  # temp threshold, below assume ppt is snow
):
    '''
    For every netatmo precipitation station,
     find nearest netatmo temperature station, intersect
     both stations and remove all days where temperature < 1°C

    '''
    # read netatmo ppt df
    in_netatmo_ppt_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                         index_col=0, sep=';',
                                         parse_dates=True,
                                         infer_datetime_format=True,
                                         engine='c')
    stns_ppt_ids = in_netatmo_ppt_stns_df.columns

    # read netatmo temp df
    in_netatmo_temperature_stns_df = pd.read_csv(netatmo_temperature_df_file,
                                                 index_col=0, sep=';',
                                                 parse_dates=True,
                                                 infer_datetime_format=True,
                                                 engine='c')

    # read netatmo distance matrix ppt-temp
    in_df_distance_netatmo_netatmo = pd.read_csv(
        distance_matrix_netatmo_ppt_netatmo_temp, sep=';', index_col=0)

    # read dwd ppt hdf5
    HDF52 = HDF5(infile=path_to_dwd_data)

    # read distance matrix dwd-netamot ppt
    in_df_distance_netatmo_dwd = pd.read_csv(
        distance_matrix_netatmo_ppt_dwd_ppt, sep=';', index_col=0)

    # read netatmo ppt coords df (for plotting)
    in_netatmo_df_coords = pd.read_csv(netatmo_ppt_coords_df, sep=';',
                                       index_col=0, engine='c')

    # create df to append results of comparing two stns
    df_results = pd.DataFrame(index=stns_ppt_ids)

    for ppt_stn_id in stns_ppt_ids:
        # iterating through netatmo ppt stations
        print('First Ppt Stn Id is', ppt_stn_id)

        # orig stn name, for locating coordinates, appending to df_results
        ppt_stn_id_name_orig = ppt_stn_id.replace('_', ':')
        try:
            # read first netatmo station
            netatmo_ppt_stn1 = in_netatmo_ppt_stns_df.loc[:, ppt_stn_id]
            netatmo_ppt_stn1.dropna(axis=0, inplace=True)
            netatmo_ppt_stn1 = netatmo_ppt_stn1[netatmo_ppt_stn1 < max_ppt_thr]

            # Flag if consider Temperature Threshold
            # remove all days with Temp < temp_thr (1°C)
            if use_temp_thr:

                # get distance netatmo ppt stn to all temp stns, get min
                distances_netatmo_to_stn1 = in_df_distance_netatmo_netatmo.loc[
                    ppt_stn_id, :]
                sorted_distances_ppt_temp = distances_netatmo_to_stn1.sort_values(
                    ascending=True)
                min_dist_ppt_temp = sorted_distances_ppt_temp.values[0]

                # check if temp station is close enough (5Km)
                if min_dist_ppt_temp <= 5000:
                    print('Seperating distance is ', min_dist_ppt_temp, ' m')

                    # if station is near, read df_temp
                    netatmo_temp_stn = sorted_distances_ppt_temp.index[0]
                    df_temp = in_netatmo_temperature_stns_df.loc[
                        :, netatmo_temp_stn]
                    df_temp.dropna(axis=0, inplace=True)

                    print(' Temp Stn Id is', netatmo_temp_stn)
                    # intersect netatmo ppt and netatmo temp data
                    idx_cmn = netatmo_ppt_stn1.index.intersection(
                        df_temp.index)
                    df_netatmo_ppt_cmn = netatmo_ppt_stn1.loc[idx_cmn]
                    df_temp_cmn = df_temp.loc[idx_cmn]

                    # if enough data available, remove all temp vals < thr
                    if (df_netatmo_ppt_cmn.values.shape[0] > 0 and
                            df_temp_cmn.values.shape[0] > 0):
                        try:
                            df_temp_abv_thr = df_temp_cmn[
                                df_temp_cmn >= temp_thr]
                            df_ppt_abv_temp_thr = df_netatmo_ppt_cmn.loc[
                                df_temp_abv_thr.index]
                        except Exception as msg:
                            print('error while selecting data abv thr', msg)
                            continue
                    else:
                        print('station is near but no common data')
                        continue
                else:
                    print('no nearby Temp station')
                    continue
            # if use temp thr and enough ppt data abv temp_thr exist
            if use_temp_thr and isinstance(df_ppt_abv_temp_thr, pd.Series):
                netatmo_ppt_stn1 = df_ppt_abv_temp_thr

            # find distance to all dwd stations, sort them, select minimum
            distances_dwd_to_stn1 = in_df_distance_netatmo_dwd.loc[
                ppt_stn_id, :]
            sorted_distances_ppt_dwd = distances_dwd_to_stn1.sort_values(
                ascending=True)
            min_dist_ppt_dwd = sorted_distances_ppt_dwd.values[0]

            if min_dist_ppt_dwd <= 5000:
                # check if dwd station is near, select and read dwd stn
                stn_2_dwd = sorted_distances_ppt_dwd.index[0]

                df_dwd = HDF52.get_pandas_dataframe(ids=[stn_2_dwd])
                df_dwd = df_dwd[df_dwd < max_ppt_thr]
                print('Second DWD Stn Id is', stn_2_dwd,
                      'distance is ', min_dist_ppt_dwd)

                # intersect dwd and netatmo ppt data
                df_netatmo_cmn, df_dwd_cmn = resample_intersect_2_dfs(
                    netatmo_ppt_stn1, df_dwd, temp_freq_resample)

                if (df_netatmo_cmn.values.shape[0] > 0 and
                        df_dwd_cmn.values.shape[0] > 0):
                    # change everything to dataframes
                    df_netatmo_cmn = pd.DataFrame(
                        data=df_netatmo_cmn.values,
                        index=df_netatmo_cmn.index,
                        columns=[ppt_stn_id])

                    df_dwd_cmn = pd.DataFrame(
                        data=df_dwd_cmn.values,
                        index=df_dwd_cmn.index,
                        columns=[stn_2_dwd])

                    # calculate prob(ppt > 1 mm) for netatmo and dwd
                    p1_netatmo = 1 - calculate_probab_ppt_below_thr(
                        df_netatmo_cmn.values, ppt_min_thr)
                    p1_dwd = 1 - calculate_probab_ppt_below_thr(
                        df_dwd_cmn.values, ppt_min_thr)
                    # compare p1 for both stns
                    (compare_p1_str,
                        colorp1) = compare_p1_dwd_p1_netatmo(
                        p1_dwd,  p1_netatmo)
                    # compare the mean of ppt for netatmo and dwd
                    (compare_mean_str,
                        colormean) = compare_p1_dwd_p1_netatmo(
                        df_dwd_cmn.values.mean(),
                        df_netatmo_cmn.values.mean())

                    # get coordinates of netatmo station
                    lon_stn_netatmo = in_netatmo_df_coords.loc[
                        ppt_stn_id_name_orig, x_col_name]
                    lat_stn_netatmo = in_netatmo_df_coords.loc[
                        ppt_stn_id_name_orig, y_col_name]

                    # append the result to df_result, for each stn
                    df_results.loc[ppt_stn_id,
                                   'lon'] = lon_stn_netatmo
                    df_results.loc[ppt_stn_id,
                                   'lat'] = lat_stn_netatmo
                    df_results.loc[ppt_stn_id,
                                   'p1'] = compare_p1_str
                    df_results.loc[ppt_stn_id,
                                   'colorp1'] = colorp1

                    df_results.loc[
                        ppt_stn_id,
                        'mean_ppt'] = compare_mean_str

                    df_results.loc[
                        ppt_stn_id,
                        'color_mean_ppt'] = colormean
                    print('ADDED DATA TO DF RESULTS')
                else:
                    print('DWD Station is near but not enough data')

            else:
                print('DWD station is not near')

        except Exception as msg:
            print('error while finding neighbours ')
            print(msg)
    # drop all netatmo stns with no data
    df_results.dropna(axis=0, how='any', inplace=True)
    return df_results


#==============================================================================
#
#==============================================================================


def plt_on_map_comparing_p1_ppt_mean_netatmo_dwd(df_results,
                                                 shp_de_file,
                                                 temp_freq,
                                                 use_temp_thr,
                                                 out_dir):
    if use_temp_thr:
        title_add = 'With_Temperatre_threshold'
    else:
        title_add = '_'

    print('plotting comparing p1')
    plt.ioff()
    fig = plt.figure(figsize=(15, 15), dpi=150)

    ax = fig.add_subplot(111)

    shp_de = shapefile.Reader(shp_de_file)
    # read and plot shapefile (BW or Germany) should be lon lat
    for shape_ in shp_de.shapeRecords():
        lon = [i[0] for i in shape_.shape.points[:][::-1]]
        lat = [i[1] for i in shape_.shape.points[:][::-1]]

        ax.scatter(lon, lat, marker='.', c='lightgrey',
                   alpha=0.25, s=2)
    # plot the stations in shapefile, look at the results of p1 comparasion
    for i in range(df_results.shape[0]):
        ax.scatter(df_results.lon.values[i],
                   df_results.lat.values[i],
                   marker=df_results.p1.values[i],
                   c=df_results.colorp1.values[i],
                   alpha=1,
                   s=15,
                   label='P1')

    ax.set_title('Difference in Probability P1 (Ppt>1mm) Netatmo and DWD %s %s'
                 % (temp_freq, title_add))
    plt.grid(alpha=0.5)
    plt.axis('equal')
    plt.savefig(
        os.path.join(
            out_dir,
            '%s_%s_differece_in_p1_netatmo_ppt_dwd_station.png'
            % (title_add, temp_freq)),

        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)
    plt.close()
    #==========================================================================
    #
    #==========================================================================
    print('plotting comparing mean')
    plt.ioff()
    fig = plt.figure(figsize=(15, 15), dpi=150)

    ax = fig.add_subplot(111)

    shp_de = shapefile.Reader(shp_de_file)

    for shape_ in shp_de.shapeRecords():
        lon = [i[0] for i in shape_.shape.points[:][::-1]]
        lat = [i[1] for i in shape_.shape.points[:][::-1]]

        ax.scatter(lon, lat, marker='.', c='lightgrey',
                   alpha=0.25, s=2)

    for i in range(df_results.shape[0]):

        ax.scatter(df_results.lon.values[i],
                   df_results.lat.values[i],
                   marker=df_results.loc[:, 'mean_ppt'].values[i],
                   c=df_results.color_mean_ppt.values[i],
                   s=15,
                   alpha=1,
                   label='Mean')
    plt.grid(alpha=0.5)
    ax.set_title('Difference in Average Rainfall values Netatmo and DWD %s %s'
                 % (temp_freq, title_add))
    plt.axis('equal')
    plt.savefig(
        os.path.join(
            out_dir,
            '%s_%s_difference_in_mean_station_ppt_dwd_station.png'
            % (title_add, temp_freq)),
        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)
    plt.close()

    return df_results

#==============================================================================
#
#==============================================================================


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program
    for temp_freq in aggregation_frequencies:
        print('Time aggregation is', temp_freq)
        df_results = select_netatmo_ppt_abv_netatmo_temp_thr(
            netatmo_ppt_df_file=path_to_ppt_netatmo_data,
            netatmo_temperature_df_file=path_to_temp_netatmo_data,
            path_to_dwd_data=path_to_ppt_hdf_data,
            netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
            distance_matrix_netatmo_ppt_netatmo_temp=distance_matrix_df_file_ppt_temp,
            distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
            temp_freq_resample=temp_freq,
            use_temp_thr=use_temp_thr,
            temp_thr=1)

        plt_on_map_comparing_p1_ppt_mean_netatmo_dwd(df_results,
                                                     path_to_shpfile,
                                                     temp_freq,
                                                     use_temp_thr,
                                                     out_save_dir_orig)
    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
