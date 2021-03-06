# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Read DWD Temperature data in BW and remove all Netatmo station
where at nearby DWD station the temperature was below 5degrees

Created on: 2019-09-12

"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# ===========================================================================
import os
import timeit
import time

import pandas as pd
import numpy as np

from scipy.spatial import cKDTree

from _00_additional_functions import resample_Humidity_Df
#==============================================================================
#
#==============================================================================
# for RH
# path_to_DWD_temp_data = (
#     r'x:\exchange\ElHachem\DWD_Reutlingen_Temperature\dwd_reutlingen_temp_2014_2019.csv')

path_to_DWD_temp_data = (
    r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\temp_dwd_2014_2019_60min.csv")

# path_to_DWD_ppt_data = (
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
#     r"\all_dwd_hourly_ppt_data_combined_2015_2019_.csv")

path_to_Netatmo_ppt_data = (
    r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\ppt_all_netatmo_rh_hourly.csv")

path_to_Netatmo_ppt_coords = (
    r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz"
    r"\netatmo_Rheinland-Pfalz_1hour_coords_utm32.csv")

path_to_DWD_temp_coords = (r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz"
                           r"\dwd_temp_coords_in_around_RH_utm32.csv")
# ppt_dwd_2014_2019_60min
# path_to_distance_mtx_Netatmo_ppt_DWD_temp = (
#     r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\distance_mtx_in_m_Netatmo_DWD_temp.csv")

# path_to_distance_mtx_Netatmo_ppt_DWD_temp = (
#     r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\distance_mtx_in_m_Netatmo_DWD_temp.csv")
#==============================================================================
#
#==============================================================================
# for BW
# path_to_DWD_temp_data = (
#     r'x:\exchange\ElHachem\DWD_Reutlingen_Temperature\dwd_reutlingen_temp_2014_2019.csv')

# path_to_DWD_temp_data = (
#     r"/run/media/abbas/EL Hachem 2019/DWD_download_Temperature_data"
#     "/all_dwd_daily_temp_data_combined_2014_2019.csv")
#
# path_to_DWD_ppt_data = (
#     r"/run/media/abbas/EL Hachem 2019/home_office/NetAtmo_BW"
#     r"/all_dwd_hourly_ppt_data_combined_2015_2019_.csv")
#
# # path_to_Netatmo_ppt_data = (
#    r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\ppt_all_netatmo_rh_hourly.csv")
# # ppt_dwd_2014_2019_60min
# # path_to_distance_mtx_Netatmo_ppt_DWD_temp = (
# #     r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\distance_mtx_in_m_Netatmo_DWD_temp.csv")
#
# # path to Netatmo ppt coords
# path_to_Netatmo_ppt_coords = (
#     r"/run/media/abbas/EL Hachem 2019/home_office/NetAtmo_BW"
#     r"/netatmo_bw_1hour_coords_utm32.csv")

# # dwd temp coords
# path_to_DWD_temp_coords = (
#     r"/run/media/abbas/EL Hachem 2019/DWD_download_Temperature_data"
#     r"/Dwd_temperature_stations_coords_in_BW_utm32.csv")
#
# path_to_DWD_ppt_coords = (
#     r"/run/media/abbas/EL Hachem 2019/home_office/NetAtmo_BW"
#     r"/station_coordinates_names_hourly_only_in_BW_utm32.csv")

freezing_temp = 5  # below 5 degrees it is snow

#==============================================================================
#
#==============================================================================


print('**** Started on %s ****\n' % time.asctime())
START = timeit.default_timer()  # to get the runtime of the program

print('Reading data\n')

# def temporal resolution of Netatmo
temp_res = '60min'

# read DWD Temp data
in_df_dwd_temp_data = pd.read_csv(
    path_to_DWD_temp_data,
    engine='c',
    sep=';', index_col=0)

# original
in_df_dwd_temp_data.index = pd.to_datetime(in_df_dwd_temp_data.index,
                                           format='%Y-%m-%d')
# in_df_dwd_temp_data.index = pd.to_datetime(in_df_dwd_temp_data.index,
#                                            format='%d-%m-%Y')

in_df_dwd_temp_data_daily = resample_Humidity_Df(in_df_dwd_temp_data,
                                                 temp_freq='D')
# read Netatmo ppt data
in_df_netatmo_ppt_data = pd.read_csv(
    path_to_Netatmo_ppt_data,  # path_to_Netatmo_ppt_data,
    sep=';', index_col=0,
    engine='c')
in_df_netatmo_ppt_data.index = pd.to_datetime(in_df_netatmo_ppt_data.index,
                                              format='%Y-%m-%d %H:%M:%S')

# read Netatmo ppt coords
in_netatmo_df_coords_utm32 = pd.read_csv(
    path_to_Netatmo_ppt_coords,  # path_to_Netatmo_ppt_coords
    sep=';', index_col=0,
    engine='c')

#in_netatmo_df_coords_utm32.index = in_df_netatmo_ppt_data.columns
# # read DWD ppt data
# in_df_dwd_ppt_data = pd.read_csv(
#     path_to_DWD_ppt_data,
#     sep=';', index_col=0,
#     engine='c')
# in_df_dwd_ppt_data.index = pd.to_datetime(in_df_dwd_ppt_data.index,
#                                               format='%Y-%m-%d %H:%M:%S')

unique_ppt_stns = in_netatmo_df_coords_utm32.index.intersection(
    in_df_netatmo_ppt_data.columns)

# read DWD temp coords
in_df_dwd_temp_coords = pd.read_csv(
    path_to_DWD_temp_coords,
    sep=';', index_col=0,
    engine='c')
# in_df_dwd_temp_coords.index = ['0' * (5 - len(str(ix))) + str(ix)
#                                for ix in in_df_dwd_temp_coords.index]


in_df_dwd_temp_coords = in_df_dwd_temp_coords.loc[
    in_df_dwd_temp_coords.index.intersection(
        in_df_dwd_temp_data.columns), :]

dwd_coords_xy = [(x, y) for x, y in zip(
    in_df_dwd_temp_coords.loc[:, 'X'].values,
    in_df_dwd_temp_coords.loc[:, 'Y'].values)]

# create a tree from coordinates
dwd_points_tree = cKDTree(dwd_coords_xy)
dwd_stns_ids = in_df_dwd_temp_coords.index
# read distance matrix dwd-netamot ppt
# in_df_distance_netatmo_dwd = pd.read_csv(
#     path_to_distance_mtx_Netatmo_ppt_DWD_temp,
#     sep=';', index_col=0)

# keep only stations in BW with data from (2014-2019)
# in_df_distance_netatmo_dwd_bw = in_df_distance_netatmo_dwd.loc[
#     in_df_netatmo_ppt_data.columns, in_df_dwd_temp_data_daily.columns]

# in_df_distance_netatmo_dwd_bw = in_df_distance_netatmo_dwd

netatmo_stn_df_no_freezing_vals = in_df_netatmo_ppt_data.loc[
    :, unique_ppt_stns]

print('Done getting all data\n')

for netatmo_stn in unique_ppt_stns:

    print('Netatmo station is', netatmo_stn)
    netatmo_stn_df = in_df_netatmo_ppt_data.loc[:, netatmo_stn].dropna()

    if netatmo_stn_df.values.size > 1:
        # find distance to all dwd stations, sort them, select minimum
        (xnetatmo, ynetamto) = (
            in_netatmo_df_coords_utm32.loc[netatmo_stn, 'X'],
            in_netatmo_df_coords_utm32.loc[netatmo_stn, 'Y'])
        # This finds the index of all points within

        distances, indices = dwd_points_tree.query(
            np.array([xnetatmo, ynetamto]),
            k=2)
    #             coords_nearest_nbr = dwd_coords_xy[indices[neighbor_to_chose]]
        dwd_stn_id = dwd_stns_ids[indices[0]]
        min_distance = np.round(distances[0], 2)

        print('Distance to DWD Temp station', np.round(min_distance, 2))
        dwd_temp_data = in_df_dwd_temp_data_daily.loc[:, dwd_stn_id].dropna()
        #.drop_duplicates()
    #     dwd_temp_data = in_df_dwd_temp_data.drop_duplicates()

        dwd_temp_data_freezing = dwd_temp_data[
            dwd_temp_data < freezing_temp].dropna()

        # freezing days and one day before
        dwd_freezing_days = dwd_temp_data_freezing.index

        print('DWD days with freezing temperatures in Netatmo Ppt',
              dwd_freezing_days.shape[0])

        # generate time period based on netatmo freq for freezing days
        intersect_ix_temp_res = []
        for freez_ix in dwd_freezing_days:

            time_range_freez = pd.date_range(start=freez_ix - pd.Timedelta(days=1),
                                             end=freez_ix +
                                             pd.Timedelta(days=1),
                                             freq=temp_res)
            intersect_ix_temp_res.append(time_range_freez)

        # make it a datetime index for intersection with netatmo station
        ix_unique_freez = pd.DatetimeIndex(
            np.unique([ix for ix_ser in intersect_ix_temp_res
                       for ix in ix_ser]))

        # find these days in netatmo ppt data
        # replace these days with nans
        intersect_ix = ix_unique_freez.intersection(netatmo_stn_df.index)

        print('%d / %d events were flagged to nans in Netatmo station'
              % (intersect_ix.size, netatmo_stn_df.size))

        netatmo_stn_df_no_freezing_vals.loc[
            intersect_ix, netatmo_stn_df.name] = np.nan

#==============================================================================
#
#==============================================================================

print('\nSaving new df')

netatmo_stn_df_no_freezing_vals.to_csv(
    r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz'
    r'\ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.csv',
    sep=';', float_format='%0.2f')
# \ppt_dwd_2014_2019_60min_no_freezing_5deg.csv'
netatmo_stn_df_no_freezing_vals.reset_index(inplace=True)
netatmo_stn_df_no_freezing_vals.rename(columns={'index': 'Time'}, inplace=True)
netatmo_stn_df_no_freezing_vals.to_feather(
    os.path.join(r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz'
                 r'\ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.fk'))

# netatmo_stn_df_no_freezing_vals.to_csv(
#     r"/run/media/abbas/EL Hachem 2019/home_office/NetAtmo_BW"
#     r"/all_dwd_hourly_ppt_data_combined_2015_2019_no_freezing_5deg_.csv",
#     sep=';', float_format='%0.2f')
# # \ppt_dwd_2014_2019_60min_no_freezing_5deg.csv'
# netatmo_stn_df_no_freezing_vals.reset_index(inplace=True)
# netatmo_stn_df_no_freezing_vals.rename(columns={'index': 'Time'}, inplace=True)
# netatmo_stn_df_no_freezing_vals.to_feather(
#     os.path.join(r"/run/media/abbas/EL Hachem 2019/home_office/NetAtmo_BW"
#                  r"/all_dwd_hourly_ppt_data_combined_2015_2019_no_freezing_5deg_.fk"))
STOP = timeit.default_timer()  # Ending time
print(('\n****Done with everything on %s.\nTotal run time was'
       ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
