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

from _00_additional_functions import (resample_intersect_2_dfs,
                                      select_convective_season,
                                      select_df_within_period,
                                      get_cdf_part_abv_thr,

                                      convert_coords_fr_wgs84_to_utm32_,
                                      )

from _10_aggregate_plot_compare_2_DWD_stns import (get_dwd_stns_coords,
                                                   get_nearest_dwd_station)

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
#==============================================================================
#
#==============================================================================

# for getting station names
path_to_ppt_netatmo_data_csv = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.csv')
assert os.path.exists(path_to_ppt_netatmo_data_csv), 'wrong NETATMO Ppt file'

# for reading ppt data station by station
# HOURLY DATA
path_to_ppt_netatmo_data_feather = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_new_no_freezing_2.fk')
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
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r"\NetAtmo_BW\all_dwd_hourly_ppt_data_combined_2014_2019_.fk")
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'

# 10Min DATA
# path_to_ppt_dwd_data = (
#     r"F:\download_DWD_data_recent\all_data_10min_ppt_data_combined_2014_2019.fk")
# assert os.path.exists(path_to_ppt_dwd_data), 'wrong 10 min DWD Ppt file'

distance_matrix_netatmo_dwd_df_file = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD.csv')
assert os.path.exists(
    distance_matrix_netatmo_dwd_df_file), 'wrong Distance MTX  file'

path_to_dwd_coords_df_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\station_coordinates_names_hourly_only_in_BW_utm32.csv")
assert os.path.exists(path_to_dwd_coords_df_file), 'wrong DWD coords file'

path_to_netatmo_coords_df_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\rain_bw_1hour"
    r"\netatmo_bw_1hour_coords.csv")
assert os.path.exists(path_to_netatmo_coords_df_file), 'wrong DWD coords file'

path_to_netatmo_gd_stns_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\plots_NetAtmo_ppt_DWD_ppt_correlation_"
    r"\keep_stns_all_neighbor_99_0_per_60min_s0.csv")
#assert os.path.exists(path_to_netatmo_gd_stns_file), 'wrong netatmo stns file'

path_to_shpfile = (r'F:\data_from_exchange\Netatmo'
                   r'\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp')

# assert os.path.exists(path_to_shpfile), 'wrong shapefile path'

out_save_dir_orig = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                     r'\plots_NetAtmo_ppt_DWD_ppt_correlation_')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

#==============================================================================
# get dwd coordinates
#==============================================================================
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

x_col_name_dwd = 'X'
y_col_name_dwd = 'Y'


in_dwd_coords_df, _, _, _ = get_dwd_stns_coords(
    path_to_dwd_coords_df_file, x_col_name_dwd, y_col_name_dwd, index_col=0,
    sep_type=';')

stndwd_ix = ['0' * (5 - len(str(stn_id))) + str(stn_id)
             if len(str(stn_id)) < 5 else str(stn_id)
             for stn_id in in_dwd_coords_df.index]

in_dwd_coords_df.index = stndwd_ix
#==============================================================================
#
#==============================================================================
# used in Netatmo coords df
x_col_name = 'lon'
y_col_name = 'lat'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 3 * 1e4  # 5000  # m

# threshold for max ppt value per hour
max_ppt_thr = 200.  # ppt above this value are not considered

# only highest x% of the values are selected
lower_percentile_val_lst = [97., 98., 99., 99.5]  # [80, 85, 90, 95, 99]


# ['10min', '60min', '120min', '480min', '720min', '1440min']
aggregation_frequencies = ['60min', '120min', '180min', '360min',
                           '720min', '1440min']

# temporal aggregation of df

# [0, 1, 2, 3, 4]  # refers to DWD neighbot (0=first)
neighbors_to_chose_lst = [0, 1, 2]  # 4, 5, 6, 7
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
end_date = '2019-09-30 00:00:00'

#==============================================================================
#
#==============================================================================


# @profile
def compare_netatmo_dwd_p1_or_p5_or_mean_ppt_or_correlations(
        path_netatmo_ppt_df_feather,  # path to df of all netatmo ppt stations
        pth_to_netatmo_cols_df_csv,  # path to csv file, get all columns
        path_to_dwd_data,  # path to dwd ppt hdf5 data
        path_netatmo_gd_stns,  # path to netatmo stns with high rank correls
        netatmo_ppt_coords_df,  # path to netatmo ppt coords df
        neighbor_to_chose,  # which DWD station neighbor to chose
        distance_matrix_netatmo_ppt_dwd_ppt,  # distance all netatmo-dwd stns
        min_dist_thr_ppt,  # distance threshold when selecting dwd neigbours
        temp_freq_resample,  # temp freq to resample dfs
        val_thr_percent,  # value in percentage, select all values above it
        min_req_ppt_vals  # threshold, minimum ppt values per station
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

    # get all station names for netatmo
#     stns_ppt_ids = pd.read_csv(
#         pth_to_netatmo_cols_df_csv, nrows=0, sep=';', engine='c',
#         memory_map=True).columns.tolist()
#     try:
#         stns_ppt_ids = list(filter(lambda x: x != 'Unnamed: 0', stns_ppt_ids))
#     except Exception as msg:
#         print(msg)

    # read distance matrix dwd-netamot ppt
    in_df_distance_netatmo_dwd = pd.read_csv(
        distance_matrix_netatmo_ppt_dwd_ppt, sep=';', index_col=0)

    # read netatmo ppt coords df (for plotting)
    in_netatmo_df_coords = pd.read_csv(netatmo_ppt_coords_df, sep=';',
                                       index_col=0, engine='c')

    # get all station names for netatmo
    stns_ppt_ids = [stn_id for stn_id in in_df_distance_netatmo_dwd.index]

    # read netatmo good stns df

    in_df_stns = pd.read_csv(path_netatmo_gd_stns, index_col=0,
                             sep=';')
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
            netatmo_ppt_stn1 = pd.read_feather(path_netatmo_ppt_df_feather,
                                               columns=['Time', ppt_stn_id],
                                               use_threads=True)
            netatmo_ppt_stn1.set_index('Time', inplace=True)
            netatmo_ppt_stn1.index = pd.to_datetime(
                netatmo_ppt_stn1.index, format=date_fmt)

#             netatmo_ppt_stn1 = in_netatmo_ppt_stns_df.loc[:, ppt_stn_id]
            netatmo_ppt_stn1.dropna(axis=0, inplace=True)
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
            distances_dwd_to_stn1 = in_df_distance_netatmo_dwd.loc[
                ppt_stn_id, :]
            sorted_distances_ppt_dwd = distances_dwd_to_stn1.sort_values(
                ascending=True)

            # select only from neighbor to chose
            sorted_distances_ppt_dwd = sorted_distances_ppt_dwd.iloc[
                neighbor_to_chose:]

            # select the DWD station neighbor
            min_dist_ppt_dwd = np.round(
                sorted_distances_ppt_dwd.values[0], 2)

            if min_dist_ppt_dwd <= min_dist_thr_ppt:
                # check if dwd station is near, select and read dwd stn
                stn_2_dwd = sorted_distances_ppt_dwd.index[0]

                df_dwd = pd.read_feather(path_to_dwd_data,
                                         columns=['Time', stn_2_dwd],
                                         use_threads=True)
                df_dwd.set_index('Time', inplace=True)

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

                        # get coordinates of netatmo station for plotting
                        lon_stn_netatmo = in_netatmo_df_coords.loc[
                            ppt_stn_id_name_orig, x_col_name]
                        lat_stn_netatmo = in_netatmo_df_coords.loc[
                            ppt_stn_id_name_orig, y_col_name]

                        #======================================================
                        # look for agreements, correlation between all values
                        #======================================================

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
                            df_netatmo_cmn >= netatmo_ppt_thr_per).astype(int)
                        df_dwd_cmn_Bool = (
                            df_dwd_cmn >= dwd_ppt_thr_per).astype(int)

                        # calculate spearman correlations of booleans 1, 0

                        bool_pears_corr = np.round(
                            pears(df_dwd_cmn_Bool.values.ravel(),
                                  df_netatmo_cmn_Bool.values.ravel())[0], 2)

#                         bool_spr_corr = np.round(
#                             spr(df_dwd_cmn_Bool.values.ravel(),
#                                 df_netatmo_cmn_Bool.values.ravel())[0], 2)
                    #==========================================================
                    # DWD stations
                    #==========================================================
                    # get id, coordinates and distances of neighbor
                        _, stn_near, distance_near = get_nearest_dwd_station(
                            first_stn_id=stn_2_dwd,
                            coords_df_file=path_to_dwd_coords_df_file,
                            x_col_name=x_col_name_dwd,
                            y_col_name=y_col_name_dwd,
                            index_col=0,
                            sep_type=';',
                            neighbor_to_chose=neighbor_to_chose + 1)

                        assert stn_2_dwd != stn_near, 'wrong DWD neighbour selected'

                        if len(stn_near) < 5:
                            stn_near = '0' * \
                                (5 - len(str(stn_near))) + str(stn_near)
                        try:

                            idf2 = pd.read_feather(path_to_ppt_dwd_data,
                                                   columns=['Time', stn_near],
                                                   use_threads=True)
                            idf2.set_index('Time', inplace=True)

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

                        print('Second DWD Stn Id is', stn_near,
                              'distance is', distance_near)

                        if distance_near < min_dist_thr_ppt:
                            print('\n resampling data')

                            if (df_dwd_cmn.values.ravel().shape[0] > min_req_ppt_vals and
                                    idf2.values.ravel().shape[0] > min_req_ppt_vals):
                                try:

                                    df_common1_dwd, df_common2_dwd = resample_intersect_2_dfs(
                                        df_dwd_cmn, idf2, temp_freq)
                                except Exception as msg:
                                    raise Exception
                            print('\n done resampling data')

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

                #                 plt.close()
                                print('data are available for correlation')

                                # calculate pearson and spearman between original
                                # values

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

                                print('\n****transform values to booleans*****\n')

                                df_dwd1_cmn_Bool = (
                                    df_common1_dwd > dwd1_ppt_thr_per).astype(int)
                                df_dwd2_cmn_Bool = (
                                    df_common2_dwd > dwd2_ppt_thr_per).astype(int)

                                # calculate spearman correlations of booleans
                                # 1, 0

            #                     bool_spr_corr = np.round(
            #                         spr(df_dwd1_cmn_Bool.values.ravel(),
            # df_dwd2_cmn_Bool.values.ravel())[0], 2)
                                bool_pears_corr_dwd = np.round(
                                    pears(df_dwd1_cmn_Bool.values.ravel(),
                                          df_dwd2_cmn_Bool.values.ravel())[0], 2)

                                if bool_pears_corr >= bool_pears_corr_dwd:
                                    print('keeping Netatmo')

                                    #==========================================
                                    # append the result to df_correlations, for each stn
                                    #==========================================
                                    df_results_correlations.loc[ppt_stn_id,
                                                                'lon'] = lon_stn_netatmo
                                    df_results_correlations.loc[ppt_stn_id,
                                                                'lat'] = lat_stn_netatmo
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
                                    print('Removing Netatmo')
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
                     'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_%d_'
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

                path_to_df_correlations = os.path.join(
                    out_save_dir_orig,
                    'new_method_pearson_year_allyears_df_comparing_correlations_max_sep_dist_%d_'
                    'freq_%s_dwd_netatmo_upper_%s_percent_data_considered'
                    '_neighbor_%d_.csv'  # filtered_95
                    % (min_dist_thr_ppt, temp_freq,
                        str(lower_percentile_val).replace('.', '_'),
                       neighbor_to_chose))

                if (not os.path.exists(path_to_df_correlations)):

                    print('\n Data frames do not exist, creating them\n')

                    (df_results_correlations
                     ) = compare_netatmo_dwd_p1_or_p5_or_mean_ppt_or_correlations(
                        path_netatmo_ppt_df_feather=path_to_ppt_netatmo_data_feather,
                        pth_to_netatmo_cols_df_csv=path_to_ppt_netatmo_data_csv,
                        path_to_dwd_data=path_to_ppt_dwd_data,
                        path_netatmo_gd_stns=path_to_netatmo_gd_stns_file,
                        netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
                        neighbor_to_chose=neighbor_to_chose,
                        distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
                        min_dist_thr_ppt=min_dist_thr_ppt,
                        temp_freq_resample=temp_freq,
                        val_thr_percent=lower_percentile_val,
                        min_req_ppt_vals=min_req_ppt_vals)
                else:
                    print('\n Data frames exist, not creating them\n')
                    df_results_correlations = pd.read_csv(path_to_df_correlations,
                                                          sep=';', index_col=0)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
