# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Filter Netatmo stations based on Indicator Correlation
Purpose: Find validity of Netatmo Station for interpolation purposes

Created on: 2019-07-16

For every Netatmo precipitation station, select the period of 2015-2019,
find the nearest DWD station, intersect both stations, same time period
construct the EDF for the Netatmo and the DWD stations respectively 
based on a given percentage threshold, find the corresponding rainfall value,
and transform all values above ppt_thr to 1 and below to 0, 2 Boolean series
calculate the pearson correlation between the two boolean stations data
for the corresponding DWD station find it's DWD neighboring station
repeat same procedure now between the DWD-DWD station and find the pearson
indicator correlation between the transformed boolean data sets.

If the value of the indicator correlation between the Netamo-DWD pair is 
greater or equal to the value betweeb the DWD-DWD pait, keep Netatmo station

Repeat this procedure for all Netatmo station, or different
quantile threshold (a probabilistic threshold) and for
different neighbors and temporal resolution.

Save the result to a df

Parameters
----------

Input Files
    DWD station data
    DWD coordinates data
    Netatmo precipitation station data
    Netatmo station coordinates data
    Distance matrix between DWD and Netatmo stations
    
Returns
-------

    
Df_correlations: df containing for every Netamo station, location (lon, lat),
    seperating distance, the pearson correlations of original data and
    for boolean transfomed data compared with the nearest DWD station.
    
Plot everything in the dataframe using a different script
Especially the change of correlation with distance
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

# generic Libs
import os
import timeit
import time

# other Libs
import pandas as pd
import numpy as np

from pathlib import Path

from pandas.plotting import register_matplotlib_converters

from scipy.stats import pearsonr as pears
from scipy.spatial import cKDTree

# own functions
from _00_additional_functions import (resample_intersect_2_dfs,
                                      select_convective_season,
                                      select_df_within_period,
                                      get_cdf_part_abv_thr,
                                      get_dwd_stns_coords,
                                      get_nearest_dwd_station)


# =============================================================================

register_matplotlib_converters()

# get working directory
# main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes')
# os.chdir(main_dir)

use_reduced_sample_dwd = False
# =============================================================================
# define path of all the needed netatmo and dwd files
# =============================================================================

# Note: for data DF, pandas feather format is used, it makes everything faster
# it allows reading data, station by station and not everything at once
# which saves memory and runs faster, changing from .csv to .fk is quite easy

# for BW
'''
# netatmo hourly data
path_to_ppt_netatmo_data_feather = Path(
    main_dir /
    r'NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.fk')
# path_to_ppt_netatmo_data_feather = Path(
#     main_dir /
#     r'NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_1deg.fk')

assert os.path.exists(
    path_to_ppt_netatmo_data_feather), 'wrong NETATMO Ppt file'

# dwd hourly data
path_to_ppt_dwd_data = Path(
    main_dir /
    r"NetAtmo_BW\all_dwd_hourly_ppt_data_combined_2014_2019_.fk")
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Ppt file'

# distance matrix netatmo-dwd
distance_matrix_netatmo_dwd_df_file = Path(
    main_dir /
    r'NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD.csv')
assert os.path.exists(
    distance_matrix_netatmo_dwd_df_file), 'wrong Distance MTX  file'

# coordinates of DWD stations (index=station names, columns=[x, y])
path_to_dwd_coords_df_file = Path(
    main_dir /
    r"NetAtmo_BW\station_coordinates_names_hourly_only_in_BW_utm32.csv")
assert os.path.exists(path_to_dwd_coords_df_file), 'wrong DWD coords file'

# coordinates of DWD stations (index=station names, columns=[x, y])
path_to_netatmo_coords_df_file = Path(
    main_dir /
    r"NetAtmo_BW\rain_bw_1hour"
    r"\netatmo_bw_1hour_coords.csv")
assert os.path.exists(
    path_to_netatmo_coords_df_file), 'wrong Netatmo coords file'

# if the filtering was done before, this allows using only remaining stations
path_to_netatmo_gd_stns_file = Path(
    main_dir /
    r"plots_NetAtmo_ppt_DWD_ppt_correlation_"
    r"keep_stns_all_neighbor_99_0_per_60min_s0.csv")
#assert os.path.exists(path_to_netatmo_gd_stns_file), 'wrong netatmo stns file'


out_save_dir_orig = (main_dir /
                     r'plots_NetAtmo_ppt_DWD_ppt_correlation_')

'''

data_dir = Path(
    r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz"
    #r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz/'
)

main_dir = data_dir

# r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz"
path_to_ppt_netatmo_data_csv = (
    data_dir / r'ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.csv')
assert os.path.exists(path_to_ppt_netatmo_data_csv), 'wrong NETATMO Ppt file'


path_to_ppt_netatmo_data_feather = (
    data_dir /
    r'ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.fk')
assert os.path.exists(
    path_to_ppt_netatmo_data_feather), 'wrong NETATMO Ppt file'


# HOURLY DATA
path_to_ppt_dwd_data = (
    data_dir /
    r'ppt_dwd_2014_2019_60min_no_freezing_5deg.fk')
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'


path_to_dwd_coords_df_file = Path(
    data_dir /
    r"dwd_coords_in_around_RH_utm32.csv")
assert os.path.exists(path_to_dwd_coords_df_file), 'wrong DWD coords file'


# distance_matrix_netatmo_dwd_df_file = (
#     data_dir /
#     r'distance_mtx_in_m_Netatmo_DWD.csv')
#
# assert os.path.exists(
#     distance_matrix_netatmo_dwd_df_file), 'wrong Distance MTX  file'

# path_to_netatmo_coords_df_file = (
#     data_dir /
#     r"netatmo_Rheinland-Pfalz_1hour_coords.csv")  # wgs84

# RLP coords in utm32
path_to_dwd_coords_utm32 = (
    # r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz'
    data_dir /
    r'dwd_coords_in_around_RH_utm32.csv')

path_to_netatmo_coords_utm32 = (
    data_dir /
    r'netatmo_Rheinland-Pfalz_1hour_coords_utm32.csv')

# path_to_netatmo_gd_stns_file = (
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
#     r"\plots_NetAtmo_ppt_DWD_ppt_correlation_"
#     r"\keep_stns_all_neighbor_99_0_per_60min_s0_comb.csv")
#assert os.path.exists(path_to_netatmo_gd_stns_file), 'wrong netatmo stns file'

# path_to_shpfile = (r'F:\data_from_exchange\Netatmo'
#                    r'\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp')

# assert os.path.exists(path_to_shpfile), 'wrong shapefile path'

# out_save_dir_orig = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#                      r'\plots_NetAtmo_ppt_DWD_ppt_correlation_reduced')
out_save_dir_orig = (
    data_dir /
    r'indicator_correlation')

if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)


# =============================================================================
if use_reduced_sample_dwd:
    # path to DWD reduced station size
    distance_matrix_netatmo_dwd_df_file = Path(
        r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
        r"\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD_reduced.csv")

    path_to_dwd_coords_df_file = Path(
        main_dir /
        r"NetAtmo_BW\reduced_smaple_dwd_stns_utm32.csv")
# =============================================================================
# get dwd coordinates
# =============================================================================
# this could be used in case coordinates needs to be transformed, this was
# done before, change Netatmo and DWD coordinates from degrees to metric
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

x_col_name_dwd = 'X'
y_col_name_dwd = 'Y'

# read dwd coords dataframe
in_dwd_coords_df, _, _, _ = get_dwd_stns_coords(
    path_to_dwd_coords_df_file, x_col_name_dwd, y_col_name_dwd,
    index_col_name='ID',
    sep_type=';')

# this is done because if a station name starts with 0, pandas reads it as an
# integer and the 0 (or 0s) are removed
# stndwd_ix = ['0' * (5 - len(str(stn_id))) + str(stn_id)
#              if len(str(stn_id)) < 5 else str(stn_id)
#              for stn_id in in_dwd_coords_df.index]
# # reindex dataframe
# in_dwd_coords_df.index = stndwd_ix
# =============================================================================
#
# =============================================================================
# used in Netatmo coords df
x_col_name = 'lon'
y_col_name = 'lat'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 100 * 1e4  # in m, for ex: 30km or 50km

# threshold for max ppt value per hour
max_ppt_thr = 200.  # ppt above this value are not considered

# only highest x% of the values are selected
lower_percentile_val_lst = [99.0]  # 97., 98., 99., 99.5

# temporal frequencies on which the filtering should be done
# , '180min', '360min', '720min', '1440min']
aggregation_frequencies = ['60min']
#
# , '120min', '180min',
#                            '360min', '720min', '1440min'

# [0, 1, 2, 3, 4]  # refers to DWD neighbot (0=first)
neighbors_to_chose_lst = [0]  # , 1, 2, 3]  # 4, 5, 6, 7

# Note: used Netatmo stations, are those with more than 2months of data
# this is done before when combining dfs

# all netatmos have more than 2 month data, this an extra check
min_req_ppt_vals = 10

# this is used to keep only data where month is not in this list
# not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april
not_convective_season = []  # this means all months are used

# date format of dataframes
date_fmt = '%Y-%m-%d %H:%M:%S'

# select data only within this period (same as netatmo)
start_date = '2015-01-01 00:00:00'
end_date = '2019-12-30 00:00:00'

# =============================================================================
# Main Function
# =============================================================================


def compare_netatmo_dwd_indicator_correlations(
        path_netatmo_ppt_df_feather,  # path to df of all netatmo ppt stations
        path_to_dwd_data,  # path to dwd ppt hdf5 data
        # path_netatmo_gd_stns,  # path to netatmo stns with high rank correls
        # netatmo_ppt_coords_df,  # path to netatmo ppt coords df
        neighbor_to_chose,  # which DWD station neighbor to chose
        # distance_matrix_netatmo_ppt_dwd_ppt,  # distance all netatmo-dwd stns
        min_dist_thr_ppt,  # distance threshold when selecting dwd neigbours
        temp_freq_resample,  # temp freq to resample dfs
        val_thr_percent,  # value in percentage, select all values above it
        min_req_ppt_vals,  # threshold, minimum ppt values per station
        path_to_dwd_coords_utm32,  # path_to_dwd_coords_utm32
        path_to_netatmo_coords_utm32  # path_to_netatmo_coords_utm32

):
    '''
     Find then for the netatmo station the neighboring DWD station
     intersect both stations, for the given probabilistic percentage
     threshold find the corresponding ppt_thr from the CDF of each station
     seperatly, make all values boolean (> 1, < 0) and calculate the pearson
     rank correlation between the two stations

     Add the result to a new dataframe and return it

    '''
    print('\n######\n getting all station names, reading dfs \n#######\n')

    # read distance matrix dwd-netamot ppt
    # in_df_distance_netatmo_dwd = pd.read_csv(
    #    distance_matrix_netatmo_ppt_dwd_ppt, sep=';', index_col=0)

    # read netatmo ppt coords df (for plotting)
    # in_netatmo_df_coords = pd.read_csv(netatmo_ppt_coords_df, sep=';',
    #                                   index_col=0, engine='c')

    # get all station names for netatmo
    #stns_ppt_ids = [stn_id for stn_id in in_df_distance_netatmo_dwd.index]

    # read netatmo good stns df
    in_netatmo_df_coords_utm32 = pd.read_csv(
        path_to_netatmo_coords_utm32, sep=';',
        index_col=0, engine='c')

    in_dwd_df_coords_utm32 = pd.read_csv(
        path_to_dwd_coords_utm32, sep=';',
        index_col=0, engine='c')

    dwd_coords_xy = [(x, y) for x, y in zip(
        in_dwd_df_coords_utm32.loc[:, 'X'].values,
        in_dwd_df_coords_utm32.loc[:, 'Y'].values)]

    # create a tree from coordinates
    dwd_points_tree = cKDTree(dwd_coords_xy)
    dwd_stns_ids = in_dwd_df_coords_utm32.index

    # get all station names for netatmo
    stns_ppt_ids = [stn_id for stn_id in in_netatmo_df_coords_utm32.index]

    # Note: why is this here?
    # read netatmo good stns df
    # could be used to check how does a filtered dataset performs on other
    # temporal frequencies or other upper percentiles

    # in_df_stns = pd.read_csv(path_netatmo_gd_stns, index_col=0,
    #                         sep=';')
    # good_stns = list(in_df_stns.values.ravel())
    # stns_ppt_ids = good_stns

    print('\n######\n done reading all dfs \n#######\n')

    # create dataframe to hold results
    df_results_correlations = pd.DataFrame(index=stns_ppt_ids)

    alls_stns_len = len(stns_ppt_ids)  # to count number of stations

    # iterating through netatmo ppt stations
    for ppt_stn_id in stns_ppt_ids:

        print('\n**\n Total number of Netatmo stations is\n**\n',
              alls_stns_len)

        alls_stns_len -= 1  # reduce number of remaining stations

#         print('\n********\n First Ppt Stn Id is', ppt_stn_id)

        # orig stn name, for locating coordinates, appending to df_results
        # Note: in orig df_coords, stn_mac is with ':' not '_'
        ppt_stn_id_name_orig = ppt_stn_id.replace('_', ':')

        try:
            # =================================================================
            # read and get all pre-work for Netatmo station
            # =================================================================

            # read first netatmo station
            netatmo_ppt_stn1 = pd.read_feather(path_netatmo_ppt_df_feather,
                                               columns=['Time', ppt_stn_id],
                                               use_threads=True)
            # convert index to datetime index
            netatmo_ppt_stn1.set_index('Time', inplace=True)
            netatmo_ppt_stn1.index = pd.to_datetime(
                netatmo_ppt_stn1.index, format=date_fmt)

            # drop nan and very very high values
            netatmo_ppt_stn1.dropna(axis=0, inplace=True)
            netatmo_ppt_stn1 = netatmo_ppt_stn1[netatmo_ppt_stn1 < max_ppt_thr]

            # select only convective season
            # only works if given list is not empty
            netatmo_ppt_stn1 = select_convective_season(
                df=netatmo_ppt_stn1,
                month_lst=not_convective_season)

            # select dataframe between recent years
            netatmo_ppt_stn1 = select_df_within_period(netatmo_ppt_stn1,
                                                       start_date,
                                                       end_date)

            # get coordinates of netatmo station for saving in DF
#             lon_stn_netatmo = in_netatmo_df_coords.loc[
#                 ppt_stn_id_name_orig, x_col_name]
#             lat_stn_netatmo = in_netatmo_df_coords.loc[
#                 ppt_stn_id_name_orig, y_col_name]

            # find distance to all dwd stations, sort them, select minimum
            (xnetatmo, ynetamto) = (
                in_netatmo_df_coords_utm32.loc[ppt_stn_id, 'X'],
                in_netatmo_df_coords_utm32.loc[ppt_stn_id, 'Y'])

            distances, indices = dwd_points_tree.query(
                np.array([xnetatmo, ynetamto]),
                k=5)

            try:
                stn_near = dwd_stns_ids[indices[neighbor_to_chose]]
            except Exception:
                print('get neighbot')
            min_dist_ppt_dwd = np.round(distances[neighbor_to_chose], 2)

#             # find distance to all DWD stations, sort them, select minimum
#             distances_dwd_to_stn1 = in_df_distance_netatmo_dwd.loc[
#                 ppt_stn_id, :]
#             sorted_distances_ppt_dwd = distances_dwd_to_stn1.sort_values(
#                 ascending=True)
#
#             # select only from neighbors to chose, defined at the begining
#             sorted_distances_ppt_dwd = sorted_distances_ppt_dwd.iloc[
#                 neighbor_to_chose:]  # this will keep the neighbor only
#
#             # get the distances to the neighboring DWD station
#             min_dist_ppt_dwd = np.round(sorted_distances_ppt_dwd.values[0], 2)

            # check if dwd station is near, if yes continue
            if min_dist_ppt_dwd <= min_dist_thr_ppt:

                # =============================================================
                # read and get all pre-work for DWD neighboring station
                # =============================================================

                # this will get the ID of the DWD station
                stn_2_dwd = stn_near  # sorted_distances_ppt_dwd.index[0]

                # TODO: FIX add TIME
                try:
                    df_dwd = pd.read_feather(path_to_dwd_data,
                                             columns=['Time', stn_2_dwd],
                                             use_threads=True)
                    df_dwd.set_index('Time', inplace=True)
                except Exception as msg:
                    #print('error reading dwd', msg)
                    df_dwd = pd.read_feather(path_to_dwd_data,
                                             columns=['index', stn_2_dwd],
                                             use_threads=True)
                    df_dwd.set_index('index', inplace=True)

                df_dwd.index = pd.to_datetime(
                    df_dwd.index, format=date_fmt)

                # remove nan and very very extreme values
                df_dwd.dropna(axis=0, inplace=True)
                df_dwd = df_dwd[df_dwd < max_ppt_thr]

                # select convective season
                df_dwd = select_convective_season(
                    df=df_dwd,
                    month_lst=not_convective_season)

                # select recent years
                df_dwd = select_df_within_period(df_dwd,
                                                 start_date,
                                                 end_date)

                # =============================================================
                # check if both stations have 'enough' data in common
                # =============================================================
                if (netatmo_ppt_stn1.values.size > min_req_ppt_vals and
                        df_dwd.values.size > min_req_ppt_vals):
                    #
                    #                     print('\n********\n Second DWD Stn Id is', stn_2_dwd,
                    #                           'distance is ', min_dist_ppt_dwd)

                    # select only data within same range
                    # DWD station data should be same time period as Netatmo
                    df_dwd = select_df_within_period(
                        df_dwd,
                        netatmo_ppt_stn1.index[0],
                        netatmo_ppt_stn1.index[-1])

                    # intersect dwd and netatmo ppt data
#                     print('intersect dwd and netatmo ppt data\n')
                    df_netatmo_cmn, df_dwd_cmn = resample_intersect_2_dfs(
                        netatmo_ppt_stn1, df_dwd, temp_freq_resample)

                    # check if intersected data exist
                    # Note: it could be that no common time period is found

                    if (df_netatmo_cmn.values.ravel().shape[0] > 0 and
                            df_dwd.values.ravel().shape[0] > 0):

                        # print('\n# Stations have common data #\n')

                        # =====================================================
                        # start the calculations
                        # =====================================================

                        # change everything to dataframes with stn Id as column
                        # this is easier for calculations
                        df_netatmo_cmn = pd.DataFrame(
                            data=df_netatmo_cmn.values,
                            index=df_netatmo_cmn.index,
                            columns=[ppt_stn_id])

                        df_dwd_cmn = pd.DataFrame(
                            data=df_dwd_cmn.values,
                            index=df_dwd_cmn.index,
                            columns=[stn_2_dwd])

                        #======================================================
                        # select only upper tail of values of both dataframes
                        #======================================================
                        val_thr_float = val_thr_percent / 100
                        # this will calculate the EDF of Netatmo station
                        netatmo_cdf_x, netatmo_cdf_y = get_cdf_part_abv_thr(
                            df_netatmo_cmn.values.ravel(), -0.1)
                        # find ppt value corresponding to quantile threshold
                        netatmo_ppt_thr_per = netatmo_cdf_x[np.where(
                            netatmo_cdf_y >= val_thr_float)][0]

                        # this will calculate the EDF of DWD station
                        dwd_cdf_x, dwd_cdf_y = get_cdf_part_abv_thr(
                            df_dwd_cmn.values.ravel(), -0.1)

                        # find ppt value corresponding to quantile threshold
                        dwd_ppt_thr_per = dwd_cdf_x[np.where(
                            dwd_cdf_y >= val_thr_float)][0]

#                         print('\n****transform values to booleans*****\n')
                        # if Xi > Ppt_thr then 1 else 0
                        df_netatmo_cmn_Bool = (
                            df_netatmo_cmn > netatmo_ppt_thr_per).astype(int)
                        df_dwd_cmn_Bool = (
                            df_dwd_cmn > dwd_ppt_thr_per).astype(int)

                        # calculate pearson correlations of booleans 1, 0

                        bool_pears_corr = np.round(
                            pears(df_dwd_cmn_Bool.values.ravel(),
                                  df_netatmo_cmn_Bool.values.ravel())[0], 2)
                    #==========================================================
                    # Check neighboring DWD stations
                    #==========================================================
                    # for the DWD station, neighboring the Netatmo station
                    # get id, coordinates and distances of DWD neighbor
                        _, stn_near, distance_near = get_nearest_dwd_station(
                            first_stn_id=stn_2_dwd,
                            coords_df_file=path_to_dwd_coords_df_file,
                            x_col_name=x_col_name_dwd,
                            y_col_name=y_col_name_dwd,
                            index_col_name='ID',
                            sep_type=';',
                            neighbor_to_chose=neighbor_to_chose + 1)

                        assert stn_2_dwd != stn_near, 'wrong DWD ngbr selected'
#                         stn_near = 'P159'
#                         if len(stn_near) < 5:
#                             stn_near = '0' * \
#                                 (5 - len(str(stn_near))) + str(stn_near)
                        try:
                            # read the neighboring DWD station
                            # TODO: FIX add TIME
                            try:
                                idf2 = pd.read_feather(path_to_dwd_data,
                                                       columns=[
                                                           'Time', stn_near],
                                                       use_threads=True)
                                idf2.set_index('Time', inplace=True)
                            except Exception as msg:
                                #print('error reading dwd', msg)
                                idf2 = pd.read_feather(path_to_dwd_data,
                                                       columns=[
                                                           'index', stn_near],
                                                       use_threads=True)
                                idf2.set_index('index', inplace=True)

                            idf2.index = pd.to_datetime(
                                idf2.index, format=date_fmt)

                            idf2.dropna(axis=0, inplace=True)

                            idf2 = select_convective_season(
                                idf2, not_convective_season)

                            idf2 = select_df_within_period(idf2,
                                                           start_date,
                                                           end_date)

                        except Exception:
                            raise Exception

#                         print('Second DWD Stn Id is', stn_near,
#                               'distance is', distance_near)
                        # calculate Indicator correlation between DWD-DWD
                        if distance_near < min_dist_thr_ppt:
                            #                             print('\n resampling data')

                            if (df_dwd_cmn.values.ravel().shape[0] >
                                min_req_ppt_vals and
                                    idf2.values.ravel().shape[0] >
                                    min_req_ppt_vals):
                                try:
                                    # resample and intersect the 2 dfs
                                    # now all 3 stations have same periods
                                    (df_common1_dwd,
                                        df_common2_dwd) = resample_intersect_2_dfs(
                                        df_dwd_cmn, idf2, temp_freq)
                                except Exception as msg:
                                    raise Exception
#                             print('\n done resampling data')
                            # same as before, now between DWD-DWD
                            if (df_common1_dwd.values.size > 0 and
                                    df_common2_dwd.values.size > 0):

                                df_common1_dwd = pd.DataFrame(
                                    data=df_common1_dwd.values,
                                    index=df_common1_dwd.index,
                                    columns=[stn_2_dwd])

                                df_common2_dwd = pd.DataFrame(
                                    data=df_common2_dwd.values,
                                    index=df_common2_dwd.index,
                                    columns=[stn_near])

#                                 print('data are available for correlation')

                                #==============================================
                                # select only upper tail of values of both dataframes
                                #==============================================
                                val_thr_float = val_thr_percent / 100

                                dwd1_cdf_x, dwd1_cdf_y = get_cdf_part_abv_thr(
                                    df_common1_dwd.values, -0.1)

                                # get dwd1 ppt thr from cdf
                                dwd1_ppt_thr_per = dwd1_cdf_x[np.where(
                                    dwd1_cdf_y >= val_thr_float)][0]

                                dwd2_cdf_x, dwd2_cdf_y = get_cdf_part_abv_thr(
                                    df_common2_dwd.values, -0.1)

                                # get dwd2 ppt thr from cdf
                                dwd2_ppt_thr_per = dwd2_cdf_x[np.where(
                                    dwd2_cdf_y >= val_thr_float)][0]

#                                 print('\n****transform values to booleans*****\n')

                                df_dwd1_cmn_Bool = (
                                    df_common1_dwd > dwd1_ppt_thr_per).astype(int)
                                df_dwd2_cmn_Bool = (
                                    df_common2_dwd > dwd2_ppt_thr_per).astype(int)

                                # calculate pearson correlations of booleans
                                # 1, 0

                                bool_pears_corr_dwd = np.round(
                                    pears(df_dwd1_cmn_Bool.values.ravel(),
                                          df_dwd2_cmn_Bool.values.ravel())[0], 2)
                                # check if indicator correlation between
                                # Netatmo and DWD is higher than between
                                # DWD and DWD neighbours, if yes, keep Netatmo
                                if bool_pears_corr >= bool_pears_corr_dwd:
                                    print('keeping Netatmo')

                                    #==========================================
                                    # append the result to df_correlations
                                    #
                                    #==========================================
#                                     df_results_correlations.loc[
#                                         ppt_stn_id,
#                                         'lon'] = lon_stn_netatmo
#                                     df_results_correlations.loc[
#                                         ppt_stn_id,
#                                         'lat'] = lat_stn_netatmo
                                    df_results_correlations.loc[
                                        ppt_stn_id,
                                        'Distance to neighbor'] = min_dist_ppt_dwd

                                    df_results_correlations.loc[
                                        ppt_stn_id,
                                        'DWD neighbor ID'] = stn_2_dwd

                                    df_results_correlations.loc[
                                        ppt_stn_id,
                                        'DWD-DWD neighbor ID'] = stn_near
                                    df_results_correlations.loc[
                                        ppt_stn_id,
                                        'Distance DWD-DWD neighbor'] = distance_near
                                    df_results_correlations.loc[
                                        ppt_stn_id,
                                        'Netatmo_%s_Per_ppt_thr'
                                        % val_thr_percent] = netatmo_ppt_thr_per

                                    df_results_correlations.loc[
                                        ppt_stn_id,
                                        'DWD_%s_Per_ppt_thr'
                                        % val_thr_percent] = dwd_ppt_thr_per

                                    df_results_correlations.loc[
                                        ppt_stn_id,
                                        'Bool_Pearson_Correlation_Netatmo_DWD'] = bool_pears_corr
                                    df_results_correlations.loc[
                                        ppt_stn_id,
                                        'Bool_Pearson_Correlation_DWD_DWD'] = bool_pears_corr_dwd
                                else:
                                    pass
                                    print('Removing Netatmo')
#                         print('\n********\n ADDED DATA TO DF RESULTS')
                    else:
                        pass
                        print('After intersecting dataframes not enough data')
                else:
                    pass
                    print('DWD Station is near but not enough data')
            else:
                pass
                print('\n********\n DWD station is not near')

        except Exception as msg:
            print('error while finding neighbours ', msg)
            raise Exception
            # continue
    assert alls_stns_len == 0, 'not all stations were considered'

    df_results_correlations.dropna(how='all', inplace=True)
    df_results_correlations.to_csv(
        os.path.join(
            out_save_dir_orig,
            'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_%d_'
            'freq_%s_dwd_netatmo_upper_%s_percent_data_considered'
            '_neighbor_%d_.csv'  # filtered_95
            % (min_dist_thr_ppt, temp_freq_resample,
               str(val_thr_percent).replace('.', '_'),
               neighbor_to_chose)),
        sep=';')

    return df_results_correlations

#==============================================================================
# CALL FUNCTION HERE
#==============================================================================


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    for lower_percentile_val in lower_percentile_val_lst:
        print('\n********\n Lower_percentile_val', lower_percentile_val)

        for temp_freq in aggregation_frequencies:
            print('\n********\n Time aggregation is', temp_freq)

            for neighbor_to_chose in neighbors_to_chose_lst:
                print('\n********\n DWD Neighbor is', neighbor_to_chose)

                # call this function to get the df, one containing
                # df_correlations comparing correlations

#                 path_to_df_correlations = os.path.join(
#                     out_save_dir_orig,
#                     'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_%d_'
#                     'freq_%s_dwd_netatmo_upper_%s_percent_data_considered'
#                     '_neighbor_%d_.csv'  # filtered_95
#                     % (min_dist_thr_ppt, temp_freq,
#                         str(lower_percentile_val).replace('.', '_'),
#                        neighbor_to_chose))
                path_to_df_correlations = ''
                if (not os.path.exists(path_to_df_correlations)):

                    print('\n Data frames do not exist, creating them\n')

                    (df_results_correlations
                     ) = compare_netatmo_dwd_indicator_correlations(
                        path_netatmo_ppt_df_feather=path_to_ppt_netatmo_data_feather,
                        path_to_dwd_data=path_to_ppt_dwd_data,
                        # path_netatmo_gd_stns=path_to_netatmo_gd_stns_file,
                        # netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
                        neighbor_to_chose=neighbor_to_chose,
                        # distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
                        min_dist_thr_ppt=min_dist_thr_ppt,
                        temp_freq_resample=temp_freq,
                        val_thr_percent=lower_percentile_val,
                        min_req_ppt_vals=min_req_ppt_vals,
                        path_to_dwd_coords_utm32=path_to_dwd_coords_utm32,
                        path_to_netatmo_coords_utm32=path_to_netatmo_coords_utm32
                    )
                else:
                    print('\n Data frames exist, not creating them\n')
                    df_results_correlations = pd.read_csv(path_to_df_correlations,
                                                          sep=';', index_col=0)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
