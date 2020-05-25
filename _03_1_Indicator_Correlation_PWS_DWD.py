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
    DWD station data
    Netatmo precipitation station data
    Netatmo station coordinates data
    
    Optional:
        # Distance matrix between DWD and Netatmo stations
        # Shapefile of BW area

Returns
-------

    
Df_correlations: df containing for every Netamo station,
    the statistical difference in terms of Pearson and Spearman 
    Correlations for original data and boolean transfomed data
    compared with the nearest DWD station.

Optional:  
    # Plot and save everything on a map using shapefile boundaries
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

from _01_2_read_hdf5 import HDF5

from _00_functions import (resample_intersect_2_dfs,
                           select_convective_season,
                           select_df_within_period,
                           get_cdf_part_abv_thr)

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
# HOURLY DATA

#==============================================================================
path_to_ppt_netatmo_data_hdf5 = (
    r"X:\staff\elhachem\2020_05_20_Netatmo_CML\netatmo_Germany__60min.h5")
assert os.path.exists(path_to_ppt_netatmo_data_hdf5), 'wrong NETATMO Ppt file'


path_to_ppt_dwd_data_hdf5 = (
    r"X:\staff\elhachem\2020_05_20_Netatmo_CML\dwd_comb_60min_data.h5")
assert os.path.exists(path_to_ppt_dwd_data_hdf5), 'wrong DWD Csv Ppt file'

#==============================================================================
# # coords of stns in utm32 in DE
#==============================================================================

# COORDINATES
path_to_dwd_coords_utm32 = (
    r"X:\staff\elhachem\2020_05_20_Netatmo_CML\DWD_60min_metadata_utm32.csv")

path_to_netatmo_coords_utm32 = (
    r"X:\staff\elhachem\2020_05_20_Netatmo_CML"
    r"\netatmo_Germany_only_1hours_coords_utm32.csv")


out_save_dir_orig = (
    r'X:\staff\elhachem\2020_05_20_Netatmo_CML\indicator_correlation')

if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

#==============================================================================
#
#==============================================================================
# used in Netatmo coords df
x_col_name = 'X'
y_col_name = 'Y'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 100 * 1e4  # 5000  # m

# threshold for max ppt value per hour
max_ppt_thr = 200.  # ppt above this value are not considered

# only highest x% of the values are selected
lower_percentile_val_lst = [99.0]  # [80, 85, 90, 95, 99]


aggregation_frequencies = ['60min']

# temporal aggregation of df

# [0, 1, 2, 3, 4]  # refers to DWD neighbot (0=first)
neighbors_to_chose_lst = [0]

# minimum hourly values that should be available per station
min_req_ppt_vals = 30 * 24 * 2

# this is used to keep only data where month is not in this list
# not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april
not_convective_season = [11, 12, 1, 2, 3]  # oct till april

date_fmt = '%Y-%m-%d %H:%M:%S'

# select data only within this period (same as netatmo)
start_date = '2018-01-01 00:00:00'
end_date = '2019-12-31 23:00:00'


#==============================================================================
#
#==============================================================================


# @profile
def indicator_correlation_netatmo_dwd(
        path_netatmo_ppt_df_hdf5,  # path to df of all netatmo ppt stations
        path_to_dwd_data_hdf5,  # path to dwd ppt hdf5 data
        path_to_netatmo_coords_utm32,  # path to netatmo coords stns utm32
        path_to_dwd_coords_utm32,  # path_to_dwd coords snts utm32
        neighbor_to_chose,  # which DWD station neighbor to chose
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
    HDF5_Netatmo = HDF5(infile=path_netatmo_ppt_df_hdf5)

    HDF5_DWD = HDF5(infile=path_to_dwd_data_hdf5)

    all_netatmo_ids = HDF5_Netatmo.get_all_names()

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

    print('\n######\n done reading all dfs \n#######\n')
    # create df to append results of comparing two stns

    df_results_correlations = pd.DataFrame(index=all_netatmo_ids)

    alls_stns_len = len(all_netatmo_ids)  # good_stns

    for ppt_stn_id in all_netatmo_ids:

        print('\n********\n Total number of Netatmo stations is\n********\n',
              alls_stns_len)
        alls_stns_len -= 1

        # iterating through netatmo ppt stations

        print('\n********\n First Ppt Stn Id is', ppt_stn_id)
        #data = HDF52.get_pandas_dataframe('P03668')
        # orig stn name, for locating coordinates, appending to df_results
        # ppt_stn_id_name_orig = ppt_stn_id.replace('_', ':')
        try:
            # read first netatmo station

            try:
                netatmo_ppt_stn1_orig = HDF5_Netatmo.get_pandas_dataframe(
                    ppt_stn_id)

            except Exception as msg:
                print('error reading dwd', msg)

            netatmo_ppt_stn1_orig = netatmo_ppt_stn1_orig[
                netatmo_ppt_stn1_orig < max_ppt_thr]

            # select df with period
            netatmo_ppt_stn1_orig = select_df_within_period(
                netatmo_ppt_stn1_orig,
                start=start_date,
                end=end_date)

            # drop all index with nan values
            netatmo_ppt_stn1_orig.dropna(axis=0, inplace=True)

            if netatmo_ppt_stn1_orig.size > min_req_ppt_vals:

                # find distance to all dwd stations, sort them, select minimum
                (xnetatmo, ynetamto) = (
                    in_netatmo_df_coords_utm32.loc[ppt_stn_id, 'X'],
                    in_netatmo_df_coords_utm32.loc[ppt_stn_id, 'Y'])

                # This finds the index of all points within
                # radius
    #             idxs_neighbours = dwd_points_tree.query_ball_point(

                distances, indices = dwd_points_tree.query(
                    np.array([xnetatmo, ynetamto]),
                    k=2)
    #             coords_nearest_nbr = dwd_coords_xy[indices[neighbor_to_chose]]
                stn_near = dwd_stns_ids[indices[neighbor_to_chose]]

                min_dist_ppt_dwd = np.round(distances[neighbor_to_chose], 2)

                if min_dist_ppt_dwd <= min_dist_thr_ppt:

                    # check if dwd station is near, select and read dwd stn
                    stn_2_dwd = stn_near
                    # TODO: FIX add TIME
                    try:
                        df_dwd_orig = HDF5_DWD.get_pandas_dataframe(stn_2_dwd)
                    except Exception as msg:
                        print('error reading dwd', msg)

                    df_dwd_orig.dropna(axis=0, inplace=True)

                    df_dwd_orig = select_df_within_period(df_dwd_orig,
                                                          start=start_date,
                                                          end=end_date)

                    # intersect dwd and netatmo ppt data

                    if (netatmo_ppt_stn1_orig.values.size > 1 and
                            df_dwd_orig.values.size > 1):
                        print('\n********\n Second DWD Stn Id is', stn_2_dwd,
                              'distance is ', min_dist_ppt_dwd)
                        # select only data within same range
                        df_dwd = select_df_within_period(
                            df_dwd_orig,
                            netatmo_ppt_stn1_orig.index[0],
                            netatmo_ppt_stn1_orig.index[-1])
                        #print('intersect dwd and netatmo ppt data\n')
                        df_netatmo_cmn, df_dwd_cmn = resample_intersect_2_dfs(
                            netatmo_ppt_stn1_orig, df_dwd_orig,
                            temp_freq_resample)

                        if (df_netatmo_cmn.values.ravel().shape[0] > 0 and
                                df_dwd.values.ravel().shape[0] > 0):
                            # print('\n# Stations have more than 2month data
                            # #\n')

                            # change everything to dataframes with stn Id
                            # as column
                            df_netatmo_cmn = pd.DataFrame(
                                data=df_netatmo_cmn.values,
                                index=df_netatmo_cmn.index,
                                columns=[ppt_stn_id])

                            df_dwd_cmn = pd.DataFrame(
                                data=df_dwd_cmn.values,
                                index=df_dwd_cmn.index,
                                columns=[stn_2_dwd])

                            # select only convective season
                            df_netatmo_cmn_season = select_convective_season(
                                df=df_netatmo_cmn,
                                month_lst=not_convective_season)

                            # select convective seasn
                            df_dwd_cmn_season = select_convective_season(
                                df=df_dwd_cmn,
                                month_lst=not_convective_season)

                            assert (df_netatmo_cmn_season.size ==
                                    df_dwd_cmn_season.size)

                            # get coordinates of netatmo station for
                            # plotting
    #                                 lon_stn_netatmo = in_netatmo_df_coords.loc[
    #                                     ppt_stn_id_name_orig, x_col_name]
    #                                 lat_stn_netatmo = in_netatmo_df_coords.loc[
    #                                     ppt_stn_id_name_orig, y_col_name]

                            #==============================================
                            # look for agreements, correlation between all values
                            #==============================================

                            # calculate pearson and spearman between original
                            # values
                            orig_pears_corr = np.round(
                                pears(df_dwd_cmn_season.values.ravel(),
                                      df_netatmo_cmn_season.values.ravel())[0], 2)

                            orig_spr_corr = np.round(
                                spr(df_dwd_cmn_season.values,
                                    df_netatmo_cmn_season.values)[0], 2)

                            #==============================================
                            # select only upper tail of values of both dataframes
                            #==============================================
                            val_thr_float = val_thr_percent / 100

                            netatmo_cdf_x, netatmo_cdf_y = get_cdf_part_abv_thr(
                                df_netatmo_cmn_season.values.ravel(), -0.1)
                            # get netatmo ppt thr from cdf
                            netatmo_ppt_thr_per = netatmo_cdf_x[np.where(
                                netatmo_cdf_y >= val_thr_float)][0]

                            dwd_cdf_x, dwd_cdf_y = get_cdf_part_abv_thr(
                                df_dwd_cmn_season.values.ravel(), -0.1)

                            # get dwd ppt thr from cdf
                            dwd_ppt_thr_per = dwd_cdf_x[np.where(
                                dwd_cdf_y >= val_thr_float)][0]

                            #print('\n****transform values to booleans*****\n')

                            df_netatmo_cmn_Bool = (
                                df_netatmo_cmn_season > netatmo_ppt_thr_per
                            ).astype(int)
                            df_dwd_cmn_Bool = (
                                df_dwd_cmn_season > dwd_ppt_thr_per
                            ).astype(int)

                            # calculate spearman correlations of booleans
                            # 1, 0

                            bool_pears_corr = np.round(
                                pears(df_dwd_cmn_Bool.values.ravel(),
                                      df_netatmo_cmn_Bool.values.ravel())[0], 2)
                            print('bool_pears_corr', bool_pears_corr)
    #                         bool_spr_corr = np.round(
    #                             spr(df_dwd_cmn_Bool.values.ravel(),
    # df_netatmo_cmn_Bool.values.ravel())[0], 2)

                            #==============================================
                            # append the result to df_correlations, for each stn
                            #==============================================
#                           df_results_correlations.loc[ppt_stn_id,
#                                                        'X'] = lon_stn_netatmo
#                           df_results_correlations.loc[ppt_stn_id,
#                                                        'Y'] = lat_stn_netatmo
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

                            print('\n********\n ADDED DATA TO DF RESULTS')
                        else:
                            print(
                                'After intersecting dataframes not enough data')
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
                     'df_indic_corr_sep_dist_%d_'
                     'freq_%s_upper_%s_per'
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

                if (not os.path.exists(path_to_df_correlations)):

                    print('\n Data frames do not exist, creating them\n')

                    (df_results_correlations
                     ) = indicator_correlation_netatmo_dwd(
                        path_netatmo_ppt_df_hdf5=path_to_ppt_netatmo_data_hdf5,
                        path_to_dwd_data_hdf5=path_to_ppt_dwd_data_hdf5,
                        path_to_netatmo_coords_utm32=path_to_netatmo_coords_utm32,
                        path_to_dwd_coords_utm32=path_to_dwd_coords_utm32,
                        neighbor_to_chose=neighbor_to_chose,
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
