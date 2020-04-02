# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Calculate and plot statistical differences between neighbours
Purpose: Find validity of Netatmo Station compared to DWD station

Created on: 2019-07-16

For every Netatmo precipitation station select the
convective season period (Mai till Ocotber), find the nearest
DWD station, intersect both stations, based on a percentage threshold,
select for every station seperatly the corresponding rainfall value based
on the CDF, using the threshold, make all values above a 1 and below a 0,
making everything boolean, calculate the spearman rank correlation between
the two stations and save the result in a df, do it considering different
neighbors and percentage threhsold ( a probabilistic threshold),
this allows capturing the change of rank correlation with distance and thr 

Do it on using all data for a station

Parameters
----------

Input Files
    DWD HDF5 station data
    Netatmo precipitation station data
    Netatmo station coordinates data
    Distance matrix between DWD and Netatmo stations
    Shapefile of BW area

Returns
-------

    
Df_correlations: df containing for every Netamo station,
    the statistical difference in terms of Pearson and Spearman 
    Correlations for original data and boolean transfomed data
    compared with the nearest DWD station.
    
Plot and save everything on a map using shapefile boundaries
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================


import os
import timeit
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from pathlib import Path
from matplotlib import rc
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears
from scipy.spatial import cKDTree

from _00_additional_functions import (resample_intersect_2_dfs,
                                      select_convective_season,
                                      select_df_within_period,
                                      select_df_within_period_year_basis,
                                      get_cdf_part_abv_thr,
                                      plt_on_map_agreements,
                                      plt_correlation_with_distance
                                      )


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

use_reduced_sample_dwd = False
#==============================================================================
#
#==============================================================================

# for getting station names in BW
# path_to_ppt_netatmo_data_csv = (
#     r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#     r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.csv')
# assert os.path.exists(path_to_ppt_netatmo_data_csv), 'wrong NETATMO Ppt file'

# path_to_ppt_netatmo_data_csv = (
#     r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz'
#     r'\ppt_all_netatmo_rh_hourly_no_freezing_5deg.csv')
# assert os.path.exists(path_to_ppt_netatmo_data_csv), 'wrong NETATMO Ppt file'

path_to_ppt_netatmo_data_csv = (
    r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz/'
    r'/ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.csv')
assert os.path.exists(path_to_ppt_netatmo_data_csv), 'wrong NETATMO Ppt file'

# for reading ppt data station by station
# HOURLY DATA

# path_to_ppt_netatmo_data_feather = (
#     r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#     r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.fk')

# path_to_ppt_netatmo_data_feather = (
#     r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz'
#     r'\ppt_all_netatmo_rh_hourly_no_freezing_5deg.fk')
# assert os.path.exists(
#     path_to_ppt_netatmo_data_feather), 'wrong NETATMO Ppt file'

path_to_ppt_netatmo_data_feather = (
    r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz'
    r'/ppt_all_netatmo_rh_2014_2019_60min_no_freezing_5deg.fk')
assert os.path.exists(
    path_to_ppt_netatmo_data_feather), 'wrong NETATMO Ppt file'


# # 10Min DATA
# path_to_ppt_netatmo_data_feather = (
#     r'F:\Netatmo_5min_data'
#     r'\ppt_all_netatmo_5min_stns_combined_.fk')
# assert os.path.exists(
#     path_to_ppt_netatmo_data_feather), 'wrong NETATMO Ppt file'

# HOURLY DATA
path_to_ppt_dwd_data = (
    #r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz'
    r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz'
    r'/ppt_dwd_2014_2019_60min_no_freezing_5deg.fk')
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'

# 10Min DATA
# path_to_ppt_dwd_data = (
#     r"F:\download_DWD_data_recent\all_data_10min_ppt_data_combined_2014_2019.fk")
# assert os.path.exists(path_to_ppt_dwd_data), 'wrong 10 min DWD Ppt file'

# BW
# distance_matrix_netatmo_dwd_df_file = (
#     r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#     r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD.csv')
# assert os.path.exists(
#     distance_matrix_netatmo_dwd_df_file), 'wrong Distance MTX  file'
#
# path_to_netatmo_coords_df_file = (
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
#     r"\rain_bw_1hour"
#     r"\netatmo_bw_1hour_coords.csv")
# assert os.path.exists(path_to_netatmo_coords_df_file), 'wrong DWD coords file'

# distance_matrix_netatmo_dwd_df_file = (
#     #r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz'
#     r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz'
#     r'/distance_mtx_in_m_Netatmo_DWD.csv')
# assert os.path.exists(
#     distance_matrix_netatmo_dwd_df_file), 'wrong Distance MTX  file'

# path_to_netatmo_coords_df_file = (
#     r"X:\staff\elhachem\Data\Netatmo_data\rain_Rheinland-Pfalz_1hour"
#     "\netatmo_Rheinland-Pfalz_1hour_coords.csv")

path_to_netatmo_coords_df_file = (
    r"/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz"
    r"/netatmo_Rheinland-Pfalz_1hour_coords_wgs84.csv")

# path_to_netatmo_gd_stns_file = (
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
#     r"\plots_NetAtmo_ppt_DWD_ppt_correlation_"
#     r"\keep_stns_all_neighbor_99_0_per_60min_s0_comb.csv")
#assert os.path.exists(path_to_netatmo_gd_stns_file), 'wrong netatmo stns file'

# path_to_shpfile = (r'F:\data_from_exchange\Netatmo'
#                    r'\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp')

# assert os.path.exists(path_to_shpfile), 'wrong shapefile path'

# COORDINATES
path_to_dwd_coords_utm32 = (
    r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz'
    r'/dwd_coords_in_around_RH_utm32.csv')
 
path_to_netatmo_coords_utm32 = (
    r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz'
    r'/netatmo_Rheinland-Pfalz_1hour_coords_utm32.csv')

# out_save_dir_orig = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#                      r'\plots_NetAtmo_ppt_DWD_ppt_correlation_reduced')
out_save_dir_orig = (
    r'/run/media/abbas/EL Hachem 2019/home_office'
    r'/2020_10_03_Rheinland_Pfalz/indicator_correlation')
#     r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\indicator_correlation')

if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

#==============================================================================

# if use_reduced_sample_dwd:
#     # path to DWD reduced station size
#     distance_matrix_netatmo_dwd_df_file = Path(
#         r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
#         r"\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD_reduced.csv")

#==============================================================================
#
#==============================================================================
# used in Netatmo coords df
x_col_name = 'lon'
y_col_name = 'lat'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 100 * 1e4  # 5000  # m

# threshold for max ppt value per hour
max_ppt_thr = 200.  # ppt above this value are not considered

# only highest x% of the values are selected
lower_percentile_val_lst = [99.0]  # [80, 85, 90, 95, 99]

# ['10min', '60min', '120min', '480min', '720min', '1440min']
aggregation_frequencies = ['60min']

# temporal aggregation of df

# [0, 1, 2, 3, 4]  # refers to DWD neighbot (0=first)
neighbors_to_chose_lst = [0]
# 30 days * 24 hours * 2month
# minimum hourly values that should be available per station
min_req_ppt_vals = 10  # 30 * 24 * 2

# this is used to keep only data where month is not in this list
# not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april
not_convective_season = []  # oct till april

plot_figures = False

date_fmt = '%Y-%m-%d %H:%M:%S'

# select data only within this period (same as netatmo)
start_date = '2015-01-01 00:00:00'
end_date = '2019-12-30 00:00:00'

#==============================================================================
#
#==============================================================================


# @profile
def compare_netatmo_dwd_p1_or_p5_or_mean_ppt_or_correlations(
        path_netatmo_ppt_df_feather,  # path to df of all netatmo ppt stations
        pth_to_netatmo_cols_df_csv,  # path to csv file, get all columns
        path_to_dwd_data,  # path to dwd ppt hdf5 data
        # path_netatmo_gd_stns,  # path to netatmo stns with high rank correls
        netatmo_ppt_coords_df,  # path to netatmo ppt coords df
        neighbor_to_chose,  # which DWD station neighbor to chose
        # distance_matrix_netatmo_ppt_dwd_ppt,  # distance all netatmo-dwd stns
        min_dist_thr_ppt,  # distance threshold when selecting dwd neigbours
        temp_freq_resample,  # temp freq to resample dfs
        val_thr_percent,  # value in percentage, select all values above it
        min_req_ppt_vals,  # threshold, minimum ppt values per station
        path_to_netatmo_coords_utm32,  # path_to_netatmo_coords_utm32
        path_to_dwd_coords_utm32,  # path_to_dwd_coords_utm32
):
    '''
     Find then for the netatmo station the neighboring DWD station
     intersect both stations, for the given probabilistic percentage
     threshold find the corresponding ppt_thr from the CDF of each station
     seperatly, make all values boolean (> 1, < 0) and calculate the spearman
     rank correlation between the two stations

     Add the result to a new dataframe and return it

    # TODO: add documentation
    '''
    print('\n######\n getting all station names, reading dfs \n#######\n')


    # read netatmo ppt coords df (for plotting)
    in_netatmo_df_coords = pd.read_csv(netatmo_ppt_coords_df, sep=';',
                                       index_col=0, engine='c')

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

#     in_df_stns = pd.read_csv(path_netatmo_gd_stns, index_col=0,
#                              sep=';')
    # good_stns = list(in_df_stns.values.ravel())
    # stns_ppt_ids = good_stns

    print('\n######\n done reading all dfs \n#######\n')
    # create df to append results of comparing two stns

    df_results_correlations = pd.DataFrame(index=stns_ppt_ids)

    alls_stns_len = len(stns_ppt_ids)  # good_stns

    for ppt_stn_id in stns_ppt_ids:

        print('\n********\n Total number of Netatmo stations is\n********\n',
              alls_stns_len)
        alls_stns_len -= 1
        # iterating through netatmo ppt stations

        print('\n********\n First Ppt Stn Id is', ppt_stn_id)

        # orig stn name, for locating coordinates, appending to df_results
        ppt_stn_id_name_orig = ppt_stn_id.replace('_', ':')
        try:
            # read first netatmo station
            try:
                netatmo_ppt_stn1 = pd.read_feather(path_netatmo_ppt_df_feather,
                                               columns=['Time', ppt_stn_id],
                                               use_threads=True)
            
                netatmo_ppt_stn1.set_index('Time', inplace=True)
            except Exception as msg:
                    # print('error reading dwd', msg)
                netatmo_ppt_stn1 = pd.read_feather(path_netatmo_ppt_df_feather,
                                         columns=['index', ppt_stn_id],
                                         use_threads=True)
                netatmo_ppt_stn1.set_index('index', inplace=True)
                
            netatmo_ppt_stn1.index = pd.to_datetime(
                netatmo_ppt_stn1.index, format=date_fmt)

            # drop all index with nan values
            netatmo_ppt_stn1.dropna(axis=0, inplace=True)
            

            # check what happens to nans
#             nans_idx = []
#             for ix in netatmo_ppt_stn1.isna().index:
#                 if netatmo_ppt_stn1.isna().loc[ix].values == True:
#                     nans_idx.append(ix)                    
#             nans_idx = pd.DatetimeIndex(nans_idx)

            # count number of zero, should remain the SAME !!!!
#             netatmo_ppt_stn1_nona = netatmo_ppt_stn1.dropna(how='all', axis=0)

#             plt.scatter(netatmo_ppt_stn1_nona.index,
#                         netatmo_ppt_stn1_nona.values)
#             nbr_zeros1 = (netatmo_ppt_stn1== 0).astype(int).sum(axis=0)
#             nbr_zeros2 = (netatmo_ppt_stn1_nona == 0).astype(int).sum(axis=0)


            netatmo_ppt_stn1 = netatmo_ppt_stn1[netatmo_ppt_stn1 < max_ppt_thr]

            # select only convective season
            netatmo_ppt_stn1 = select_convective_season(
                df=netatmo_ppt_stn1,
                month_lst=not_convective_season)
            # select df recent years
            netatmo_ppt_stn1 = select_df_within_period(netatmo_ppt_stn1,
                                                       start_date,
                                                       end_date)

            
            # find distance to all dwd stations, sort them, select minimum
            (xnetatmo, ynetamto) = (
                in_netatmo_df_coords_utm32.loc[ppt_stn_id, 'X'],
                in_netatmo_df_coords_utm32.loc[ppt_stn_id, 'Y'])
            # This finds the index of all points within
            
            distances, indices = dwd_points_tree.query(
                np.array([xnetatmo, ynetamto]),
                       k=5)
            
            stn_near = dwd_stns_ids[indices[neighbor_to_chose]]
            
            # for BW
            # stn_near = '0' * (5 - len(str(stn_near))) + str(stn_near)
            
            min_dist_ppt_dwd = np.round(distances[neighbor_to_chose], 2)


            if min_dist_ppt_dwd <= min_dist_thr_ppt:
                # check if dwd station is near, select and read dwd stn
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

                df_dwd.dropna(axis=0, inplace=True)

                df_dwd = df_dwd[df_dwd < max_ppt_thr]

                # select convective season
                df_dwd = select_convective_season(
                    df=df_dwd,
                    month_lst=not_convective_season)

                df_dwd = select_df_within_period(df_dwd,
                                                 start_date,
                                                 end_date)
                # intersect dwd and netatmo ppt data

                if (netatmo_ppt_stn1.values.size > min_req_ppt_vals and
                        df_dwd.values.size > min_req_ppt_vals):
                    print('\n********\n Second DWD Stn Id is', stn_2_dwd,
                          'distance is ', min_dist_ppt_dwd)
                    # select only data within same range
                    df_dwd = select_df_within_period(df_dwd,
                                                     netatmo_ppt_stn1.index[0],
                                                     netatmo_ppt_stn1.index[-1])
                    print('intersect dwd and netatmo ppt data\n')
                    df_netatmo_cmn, df_dwd_cmn = resample_intersect_2_dfs(
                        netatmo_ppt_stn1, df_dwd, temp_freq_resample)

                    if (df_netatmo_cmn.values.ravel().shape[0] > 0 and
                            df_dwd.values.ravel().shape[0] > 0):
                        print('\n# Stations have more than 2month data #\n')

                        # change everything to dataframes with stn Id as column
                        df_netatmo_cmn = pd.DataFrame(
                            data=df_netatmo_cmn.values,
                            index=df_netatmo_cmn.index,
                            columns=[ppt_stn_id])

                        df_dwd_cmn = pd.DataFrame(
                            data=df_dwd_cmn.values,
                            index=df_dwd_cmn.index,
                            columns=[stn_2_dwd])
                                                
                        #======================================================
                        # # shift only part
                        #======================================================
                        
#                         in_df_april_mid_oct = select_df_within_period_year_basis(
#                             df_netatmo_cmn)
#     
#                         # shift by one hour forward
#                         in_df_april_mid_oct = in_df_april_mid_oct.shift(1)
#                         
#                         idx_to_keep = [
#                             ix for ix in df_netatmo_cmn.index
#                              if ix not in in_df_april_mid_oct.index]
#                         in_df_mid_oct_mars = df_netatmo_cmn.loc[idx_to_keep, :]
#                         
#                         time_arr = pd.date_range(start=df_netatmo_cmn.index[0],
#                                                 end=df_netatmo_cmn.index[-1],
#                                                 freq='H')
#                         
#                         df_shifted = pd.DataFrame(index=time_arr,
#                                                 columns=['sum_rain'])
#                         df_shifted.loc[
#                             in_df_april_mid_oct.index,
#                              'sum_rain'
#                              ] = in_df_april_mid_oct.values.ravel()
#                              
#                         df_shifted.loc[
#                             in_df_mid_oct_mars.index,
#                             'sum_rain'] = in_df_mid_oct_mars.values.ravel()
#                                                  
#                         df_shifted.dropna(inplace=True, how='all')
                        df_shifted = df_netatmo_cmn.shift(1)
                        
                        df_netatmo_cmn_shifted = df_shifted 
                        cmn_idx = df_dwd_cmn.index.intersection(
                            df_netatmo_cmn_shifted.index)
                        df_dwd_cmn = df_dwd_cmn.loc[cmn_idx, :]
                        
                        df_netatmo_cmn = df_netatmo_cmn.loc[cmn_idx, :]
                        #======================================================
                        # 
                        #======================================================
                        # get coordinates of netatmo station for plotting
                        # lon_stn_netatmo = in_netatmo_df_coords.loc[
                        #    ppt_stn_id_name_orig, x_col_name]
                        # lat_stn_netatmo = in_netatmo_df_coords.loc[
                        #    ppt_stn_id_name_orig, y_col_name]

                        #======================================================
                        # look for agreements, correlation between all values
                        #======================================================

                        # calculate pearson and spearman between original
                        # values
                        orig_pears_corr = np.round(
                            pears(df_dwd_cmn.values.ravel(),
                                  df_netatmo_cmn.values.ravel())[0], 2)

                        orig_spr_corr = np.round(
                            spr(df_dwd_cmn.values,
                                df_netatmo_cmn.values)[0], 2)

                        #======================================================
                        # select only upper tail of values of both dataframes
                        #======================================================
                        val_thr_float = val_thr_percent / 100

                        netatmo_cdf_x, netatmo_cdf_y = get_cdf_part_abv_thr(
                            df_netatmo_cmn.values.ravel(), -0.1)
                        # get netatmo ppt thr from cdf
                        netatmo_ppt_thr_per = netatmo_cdf_x[np.where(
                            netatmo_cdf_y >= val_thr_float)][0]

                        dwd_cdf_x, dwd_cdf_y = get_cdf_part_abv_thr(
                            df_dwd_cmn.values.ravel(), -0.1)

                        # get dwd ppt thr from cdf
                        dwd_ppt_thr_per = dwd_cdf_x[np.where(
                            dwd_cdf_y >= val_thr_float)][0]

                        print('\n****transform values to booleans*****\n')

                        df_netatmo_cmn_Bool = (
                            df_netatmo_cmn > netatmo_ppt_thr_per).astype(int)
                            
                        df_netatmo_cmn_Bool_shifted = (
                            df_netatmo_cmn_shifted > netatmo_ppt_thr_per).astype(int)

                        df_dwd_cmn_Bool = (
                            df_dwd_cmn > dwd_ppt_thr_per).astype(int)

                        # calculate spearman correlations of booleans 1, 0

                        bool_pears_corr = np.round(
                            pears(df_dwd_cmn_Bool.values.ravel(),
                                  df_netatmo_cmn_Bool.values.ravel())[0], 4)
                        
                        bool_pears_corr_shifted = np.round(
                            pears(df_dwd_cmn_Bool.values.ravel(),
                                  df_netatmo_cmn_Bool_shifted.values.ravel())[0], 4)
                        
#                         bool_spr_corr = np.round(
#                             spr(df_dwd_cmn_Bool.values.ravel(),
#                                 df_netatmo_cmn_Bool.values.ravel())[0], 2)

                        #======================================================
                        # append the result to df_correlations, for each stn
                        #======================================================
#                         df_results_correlations.loc[ppt_stn_id,
#                                                     'lon'] = lon_stn_netatmo
#                         df_results_correlations.loc[ppt_stn_id,
#                                                     'lat'] = lat_stn_netatmo
                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Distance to neighbor'] = min_dist_ppt_dwd

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'DWD neighbor ID'] = stn_2_dwd

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
                            'Orig_Pearson_Correlation'] = orig_pears_corr

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Orig_Spearman_Correlation'] = orig_spr_corr

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Bool_Spearman_Correlation'] = bool_pears_corr
                            
                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Bool_Spearman_Correlation_Shifted'
                            ] = bool_pears_corr_shifted

                        print('\n********\n ADDED DATA TO DF RESULTS')
                    else:
                        print('After intersecting dataframes not enough data')
                else:
                    print('DWD Station is near but not enough data')
            else:
                print('\n********\n DWD station is not near')

        except Exception as msg:
            print('error while finding neighbours ', msg)
            continue
    assert alls_stns_len == 0, 'not all stations were considered'

    df_results_correlations.dropna(how='all', inplace=True)
    df_results_correlations.to_csv(
        os.path.join(out_save_dir_orig,
                     '4pearson_year_allyears_df_comparing_correlations_max_sep_dist_%d_'
                     'freq_%s_dwd_netatmo_upper_%s_percent_data_considered'
                     '_neighbor_%d_.csv'  # filtered_95
                     % (min_dist_thr_ppt, temp_freq_resample,
                        str(val_thr_percent).replace('.', '_'),
                         neighbor_to_chose)),
        sep=';')

    return df_results_correlations

#==============================================================================
#
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
                path_to_df_correlations = ''
#                 path_to_df_correlations = os.path.join(
#                     out_save_dir_orig,
#                     'pearson_year_allyears_df_comparing_correlations_max_sep_dist_%d_'
#                     'freq_%s_dwd_netatmo_upper_%s_percent_data_considered'
#                     '_neighbor_%d_filtered_99.csv'  # filtered_95
#                     % (min_dist_thr_ppt, temp_freq,
#                         str(lower_percentile_val).replace('.', '_'),
#                        neighbor_to_chose))

                if (not os.path.exists(path_to_df_correlations)):

                    print('\n Data frames do not exist, creating them\n')

                    (df_results_correlations
                     ) = compare_netatmo_dwd_p1_or_p5_or_mean_ppt_or_correlations(
                        path_netatmo_ppt_df_feather=path_to_ppt_netatmo_data_feather,
                        pth_to_netatmo_cols_df_csv=path_to_ppt_netatmo_data_csv,
                        path_to_dwd_data=path_to_ppt_dwd_data,
                        # path_netatmo_gd_stns=path_to_netatmo_gd_stns_file,
                        netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
                        neighbor_to_chose=neighbor_to_chose,
                        # distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
                        min_dist_thr_ppt=min_dist_thr_ppt,
                        temp_freq_resample=temp_freq,
                        val_thr_percent=lower_percentile_val,
                        min_req_ppt_vals=min_req_ppt_vals,
                        path_to_netatmo_coords_utm32=path_to_netatmo_coords_utm32,
                        path_to_dwd_coords_utm32=path_to_dwd_coords_utm32)
                else:
                    print('\n Data frames exist, not creating them\n')
                    df_results_correlations = pd.read_csv(path_to_df_correlations,
                                                          sep=';', index_col=0)

#                 if plot_figures:
#                     print('\n********\n Plotting Correlation with distance')
#                     plt_correlation_with_distance(
#                         df_correlations=df_results_correlations,
#                         dist_col_to_plot='Distance to neighbor',
#                         corr_col_to_plot='Bool_Spearman_Correlation',
#                         temp_freq=temp_freq,
#                         out_dir=out_save_dir_orig,
#                         year_vals='all_years',
#                         val_thr_percent=lower_percentile_val,
#                         neighbor_nbr=neighbor_to_chose)
# 
#                     plt_correlation_with_distance(
#                         df_correlations=df_results_correlations,
#                         dist_col_to_plot='Distance to neighbor',
#                         corr_col_to_plot='Orig_Spearman_Correlation',
#                         temp_freq=temp_freq,
#                         out_dir=out_save_dir_orig,
#                         year_vals='all_years',
#                         val_thr_percent=lower_percentile_val,
#                         neighbor_nbr=neighbor_to_chose)
# 
#                     print('\n********\n Plotting Correlation maps')
#                     for col_label in df_results_correlations.columns:
#                         if ('Correlation' in col_label):
#                                 # and 'Bool_Spearman' in col_label):
#                             # plot the results of df_results_correlations
#                             plt_on_map_agreements(
#                                 df_correlations=df_results_correlations,
#                                 col_to_plot=col_label,
#                                 shp_de_file=path_to_shpfile,
#                                 temp_freq=temp_freq,
#                                 out_dir=out_save_dir_orig,
#                                 year_vals=('all_years_%d_m_distance_neighbor_%d_'
#                                            % (min_dist_thr_ppt, neighbor_to_chose)),
#                                 val_thr_percent=lower_percentile_val)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
