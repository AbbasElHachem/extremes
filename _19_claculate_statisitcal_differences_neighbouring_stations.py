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


import os
import timeit
import time
import shapefile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from adjustText import adjust_text
from pathlib import Path
from matplotlib import rc
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters
from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

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
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_ppt_Netatmo_temp.csv')
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

path_to_shpfile = (r'X:\exchange\ElHachem\Netatmo'
                   r'\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp')

assert os.path.exists(path_to_shpfile), 'wrong shapefile path'

out_save_dir_orig = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                     r'\plots_NetAtmo_ppt_NetAtmo_temperature')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

#==============================================================================
#
#==============================================================================
x_col_name = ' lon'
y_col_name = ' lat'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 5000  # m
min_dist_thr_temp = 5000  # m

# threshold for max ppt value per hour
max_ppt_thr = 100.
ppt_min_thr = 1  # used when calculating p1 = 1-p0

# aggregation_frequencies = ['60min', '120min', '180min', '240min',
#                            '360min', '720min', '1440min']
aggregation_frequencies = ['60min']

# if True remove all Ppt values where Temp < Temp_thr= 1°C
temp_thr = 1
use_temp_thr = False

# this is used to keep only data where month is not in this list
not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april

# if True, work year by year, to find tendency of a stations
year_by_year = True
#==============================================================================
#
#==============================================================================


def select_convective_season(
    df,  # df to slice
    month_lst  # list of month for convective season
):
    '''
    return dataframe without the data corresponding to the winter season
    '''
    df = df.copy()
    df_conv_season = df[~df.index.month.isin(month_lst)]

    return df_conv_season


#==============================================================================
#
#==============================================================================

def compare_p1_dwd_p1_netatmo(
    val_dwd,  # p1 or ppt_mean  for df_DWD
    val_netatmo  # p1 or ppt_mean for df_Netatmo
):
    '''
    use this function to compare two values from two different stations
    find if the values are equal +-10% , if smaller or if bigger
    return the result for being saved in a dataframe
    plot the result as scatter plots

    Return:
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


def look_agreement_df1_df2(df_dwd,  # df dwd intersected with netatmo
                           df_netatmo,  # df netatmo intersected with dwd
                           ppt_thr,  # ppt thr to check agreement
                           ):
    '''
    read two dataframes, one for dwd and one for netatmo
    replace all values above threshold with 1 and below with 0

    calculate coorelation between the two
    look for agreement or disagreement
    try it with or without considering temperature threshold
    '''
    df_dwd_copy = df_dwd.copy()
    df_netatmo_copy = df_netatmo.copy()
    # calculate pearson and spearman between original values

#     df_dwd_copy_abv_thr = df_dwd_copy[df_dwd_copy > ppt_thr]
#     df_netatmo_copy_abv_thr = df_netatmo_copy[df_netatmo_copy > ppt_thr]

    orig_pear_corr = np.round(
        pears(df_dwd_copy.values.ravel(),
              df_netatmo_copy.values.ravel())[0], 2)

    orig_spr_corr = np.round(spr(df_dwd_copy.values.ravel(),
                                 df_netatmo_copy.values.ravel())[0], 2)
    # transform values to booleans
    df_dwd_copy['Bool'] = (df_dwd_copy.values > ppt_thr).astype(int)
    df_netatmo_copy['Bool'] = (df_netatmo_copy.values > ppt_thr).astype(int)

    # calculate correlations of booleans 1, 0 (pearson and spearman)
    pearson_corr = np.round(pears(df_dwd_copy.Bool.values.ravel(),
                                  df_netatmo_copy.Bool.values.ravel())[0], 2)

    spearman_corr = np.round(spr(df_dwd_copy.Bool.values.ravel(),
                                 df_netatmo_copy.Bool.values.ravel())[0], 2)

    return orig_pear_corr, orig_spr_corr, pearson_corr, spearman_corr
#==============================================================================
#
#==============================================================================


def compare_netatmo_dwd_p1_or_p5_or_mean_ppt(
        netatmo_ppt_df_file,  # path to df of all netatmo ppt stations
        path_to_dwd_data,  # path to dwd ppt hdf5 data
        netatmo_ppt_coords_df,  # path to netatmo ppt coords df
        distance_matrix_netatmo_ppt_dwd_ppt,
        min_dist_thr_ppt,  # distance threshold when selecting dwd neigbours
        temp_freq_resample,  # temp freq to resample dfs
        df_min_ppt_thr,  # ppt_thr, select all values above thr
):
    '''
    For every netatmo precipitation station,
     find nearest netatmo temperature station, intersect
     both stations and remove all days where temperature < 1°C

     Find then for the netatmo station the neares DWD station
     intersect both stations, calculate the difference in the
     probability that P>1mm and in the mean value

     Add the result to a new dataframe and return it

    '''
    print('\n######\n reading all dfs \n#######\n')
    # read netatmo ppt df
    in_netatmo_ppt_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                         index_col=0, sep=';',
                                         parse_dates=True,
                                         infer_datetime_format=True,
                                         engine='c')
    stns_ppt_ids = in_netatmo_ppt_stns_df.columns

    # read dwd ppt hdf5
    HDF52 = HDF5(infile=path_to_dwd_data)

    # read distance matrix dwd-netamot ppt
    in_df_distance_netatmo_dwd = pd.read_csv(
        distance_matrix_netatmo_ppt_dwd_ppt, sep=';', index_col=0)

    # read netatmo ppt coords df (for plotting)
    in_netatmo_df_coords = pd.read_csv(netatmo_ppt_coords_df, sep=';',
                                       index_col=0, engine='c')
    print('\n######\n done reading all dfs \n#######\n')
    # create df to append results of comparing two stns
    df_results = pd.DataFrame(index=stns_ppt_ids)
    df_results_correlations = pd.DataFrame(index=stns_ppt_ids)

    alls_stns_len = len(stns_ppt_ids)
    for ppt_stn_id in stns_ppt_ids:
        print('\n********\n Total number of Netatmo stations is\n********\n',
              alls_stns_len)
        # iterating through netatmo ppt stations

        print('\n********\n First Ppt Stn Id is', ppt_stn_id)

        # orig stn name, for locating coordinates, appending to df_results
        ppt_stn_id_name_orig = ppt_stn_id.replace('_', ':')
        try:
            # read first netatmo station
            netatmo_ppt_stn1 = in_netatmo_ppt_stns_df.loc[:, ppt_stn_id]
            netatmo_ppt_stn1.dropna(axis=0, inplace=True)
            netatmo_ppt_stn1 = netatmo_ppt_stn1[netatmo_ppt_stn1 < max_ppt_thr]

            # select only convective season
            netatmo_ppt_stn1 = select_convective_season(
                df=netatmo_ppt_stn1,
                month_lst=not_convective_season)

            # find distance to all dwd stations, sort them, select minimum
            distances_dwd_to_stn1 = in_df_distance_netatmo_dwd.loc[
                ppt_stn_id, :]
            sorted_distances_ppt_dwd = distances_dwd_to_stn1.sort_values(
                ascending=True)
            min_dist_ppt_dwd = sorted_distances_ppt_dwd.values[0]

            if min_dist_ppt_dwd <= min_dist_thr_ppt:
                # check if dwd station is near, select and read dwd stn
                stn_2_dwd = sorted_distances_ppt_dwd.index[0]

                df_dwd = HDF52.get_pandas_dataframe(ids=[stn_2_dwd])
                df_dwd = df_dwd[df_dwd < max_ppt_thr]

                # select convective season
                df_dwd = select_convective_season(
                    df=df_dwd,
                    month_lst=not_convective_season)

                print('\n********\n Second DWD Stn Id is', stn_2_dwd,
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

                    # get coordinates of netatmo station
                    lon_stn_netatmo = in_netatmo_df_coords.loc[
                        ppt_stn_id_name_orig, x_col_name]
                    lat_stn_netatmo = in_netatmo_df_coords.loc[
                        ppt_stn_id_name_orig, y_col_name]

                    #==================================================
                    #  compare p1 and mean ppt for both stns
                    #==================================================
                    # calculate prob(ppt > 1 mm) for netatmo and dwd
                    p1_netatmo = 1 - calculate_probab_ppt_below_thr(
                        df_netatmo_cmn.values, df_min_ppt_thr)
                    p1_dwd = 1 - calculate_probab_ppt_below_thr(
                        df_dwd_cmn.values, df_min_ppt_thr)
                    # compare p1 for both stns
                    (compare_p1_str,
                        colorp1) = compare_p1_dwd_p1_netatmo(
                        p1_dwd,  p1_netatmo)
                    # compare the mean of ppt for netatmo and dwd
                    (compare_mean_str,
                        colormean) = compare_p1_dwd_p1_netatmo(
                        df_dwd_cmn.values.mean(),
                        df_netatmo_cmn.values.mean())

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
                    #==================================================
                    # look for agreements
                    #==================================================
                    # calculate correlation, look for agreements
                    (orig_pears_corr, orig_spr_coor,
                        bool_pears_corr, bool_spr_corr) = look_agreement_df1_df2(
                        df_dwd=df_dwd_cmn,
                        df_netatmo=df_netatmo_cmn,
                        ppt_thr=df_min_ppt_thr)

                    # append the result to df_result, for each stn
                    df_results_correlations.loc[ppt_stn_id,
                                                'lon'] = lon_stn_netatmo
                    df_results_correlations.loc[ppt_stn_id,
                                                'lat'] = lat_stn_netatmo
                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Orig Pearson Correlation'] = orig_pears_corr

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Orig Spearman Correlation'] = orig_spr_coor

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Bool Pearson Correlation'] = bool_pears_corr

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Bool Spearman Correlation'] = bool_spr_corr

                    print('\n********\n ADDED DATA TO DF RESULTS')

                else:
                    print('DWD Station is near but not enough data')
            else:
                print('\n********\n DWD station is not near')

        except Exception as msg:
            print('error while finding neighbours ', msg)
        alls_stns_len -= 1

    # drop all netatmo stns with no data
    df_results.dropna(axis=0, how='any', inplace=True)
    df_results_correlations.dropna(axis=0, how='any', inplace=True)

    # save both dataframes
    df_results.to_csv(
        os.path.join(
            out_save_dir_orig,
            'df_comparing_p%d_and_mean_freq_%s_dwd_netatmo_with_temp_thr_%s.csv'
            % (df_min_ppt_thr, temp_freq_resample, use_temp_thr)),
        sep=';')
    df_results_correlations.to_csv(
        os.path.join(
            out_save_dir_orig,
            'df_comparing_correlations_p%d_freq_%s_dwd_netatmo_with_temp_thr_%s.csv'
            % (df_min_ppt_thr, temp_freq_resample, use_temp_thr)),
        sep=';')
    return df_results, df_results_correlations

#==============================================================================
#
#==============================================================================


def compare_netatmo_dwd_p1_or_p5_or_mean_ppt_year_by_year(
        netatmo_ppt_df_file,  # path to df of all netatmo ppt stations
        path_to_dwd_data,  # path to dwd ppt hdf5 data
        netatmo_ppt_coords_df,  # path to netatmo ppt coords df
        distance_matrix_netatmo_ppt_dwd_ppt,
        min_dist_thr_ppt,  # distance threshold when selecting dwd neigbours
        temp_freq_resample,  # temp freq to resample dfs
        df_min_ppt_thr,  # ppt_thr, select all values above thrs
        year_val,  # calculate only for this year
):
    '''
    For every netatmo precipitation station,
     find nearest netatmo temperature station, intersect
     both stations and remove all days where temperature < 1°C

     Find then for the netatmo station the neares DWD station
     intersect both stations, calculate the difference in the
     probability that P>1mm and in the mean value

     Add the result to a new dataframe and return it

    '''
    print('\n######\n reading all dfs \n#######\n')
    # read netatmo ppt df
    in_netatmo_ppt_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                         index_col=0, sep=';',
                                         parse_dates=True,
                                         infer_datetime_format=True,
                                         engine='c')
    stns_ppt_ids = in_netatmo_ppt_stns_df.columns

    # read dwd ppt hdf5
    HDF52 = HDF5(infile=path_to_dwd_data)

    # read distance matrix dwd-netamot ppt
    in_df_distance_netatmo_dwd = pd.read_csv(
        distance_matrix_netatmo_ppt_dwd_ppt, sep=';', index_col=0)

    # read netatmo ppt coords df (for plotting)
    in_netatmo_df_coords = pd.read_csv(netatmo_ppt_coords_df, sep=';',
                                       index_col=0, engine='c')
    print('\n######\n done reading all dfs \n#######\n')
    # create df to append results of comparing two stns
    df_results = pd.DataFrame(index=stns_ppt_ids)
    df_results_correlations = pd.DataFrame(index=stns_ppt_ids)

    alls_stns_len = len(stns_ppt_ids)
    for ppt_stn_id in stns_ppt_ids:
        print('\n********\n Total number of Netatmo stations is\n********\n',
              alls_stns_len)
        # iterating through netatmo ppt stations

        print('\n********\n First Ppt Stn Id is', ppt_stn_id)

        # orig stn name, for locating coordinates, appending to df_results
        ppt_stn_id_name_orig = ppt_stn_id.replace('_', ':')
        try:
            # read first netatmo station
            netatmo_ppt_stn1 = in_netatmo_ppt_stns_df.loc[:, ppt_stn_id]
            netatmo_ppt_stn1.dropna(axis=0, inplace=True)
            netatmo_ppt_stn1 = netatmo_ppt_stn1[netatmo_ppt_stn1 < max_ppt_thr]

            # select only convective season
            netatmo_ppt_stn1 = select_convective_season(
                df=netatmo_ppt_stn1,
                month_lst=not_convective_season)

            # find distance to all dwd stations, sort them, select minimum
            distances_dwd_to_stn1 = in_df_distance_netatmo_dwd.loc[
                ppt_stn_id, :]
            sorted_distances_ppt_dwd = distances_dwd_to_stn1.sort_values(
                ascending=True)
            min_dist_ppt_dwd = sorted_distances_ppt_dwd.values[0]

            if int(year_val) in np.unique(netatmo_ppt_stn1.index.year):
                print('\n++\n Year is', year_val, '\n++\n')
                netatmo_ppt_stn1 = netatmo_ppt_stn1.loc[
                    netatmo_ppt_stn1.index.year == int(year_val)]

                if min_dist_ppt_dwd <= min_dist_thr_ppt:
                    # check if dwd station is near, select and read dwd stn
                    stn_2_dwd = sorted_distances_ppt_dwd.index[0]

                    df_dwd = HDF52.get_pandas_dataframe(ids=[stn_2_dwd])
                    df_dwd = df_dwd[df_dwd < max_ppt_thr]

                    # select convective season
                    df_dwd = select_convective_season(
                        df=df_dwd,
                        month_lst=not_convective_season)

                    print('\n********\n Second DWD Stn Id is', stn_2_dwd,
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

                        # get coordinates of netatmo station
                        lon_stn_netatmo = in_netatmo_df_coords.loc[
                            ppt_stn_id_name_orig, x_col_name]
                        lat_stn_netatmo = in_netatmo_df_coords.loc[
                            ppt_stn_id_name_orig, y_col_name]

                        #==================================================
                        #  compare p1 and mean ppt for both stns
                        #==================================================
                        # calculate prob(ppt > 1 mm) for netatmo and dwd
                        p1_netatmo = 1 - calculate_probab_ppt_below_thr(
                            df_netatmo_cmn.values, df_min_ppt_thr)
                        p1_dwd = 1 - calculate_probab_ppt_below_thr(
                            df_dwd_cmn.values, df_min_ppt_thr)
                        # compare p1 for both stns
                        (compare_p1_str,
                            colorp1) = compare_p1_dwd_p1_netatmo(
                            p1_dwd,  p1_netatmo)
                        # compare the mean of ppt for netatmo and dwd
                        (compare_mean_str,
                            colormean) = compare_p1_dwd_p1_netatmo(
                            df_dwd_cmn.values.mean(),
                            df_netatmo_cmn.values.mean())

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
                        #==================================================
                        # look for agreements
                        #==================================================
                        # calculate correlation, look for agreements
                        (orig_pears_corr,
                         orig_spr_coor,
                         bool_pears_corr,
                         bool_spr_corr) = look_agreement_df1_df2(
                            df_dwd=df_dwd_cmn,
                            df_netatmo=df_netatmo_cmn,
                            ppt_thr=df_min_ppt_thr)

                        # append the result to df_result, for each stn
                        df_results_correlations.loc[
                            ppt_stn_id, 'lon'] = lon_stn_netatmo
                        df_results_correlations.loc[
                            ppt_stn_id, 'lat'] = lat_stn_netatmo
                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Orig Pearson Correlation'] = orig_pears_corr

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Orig Spearman Correlation'] = orig_spr_coor

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Bool Pearson Correlation'] = bool_pears_corr

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Bool Spearman Correlation'] = bool_spr_corr

                        print('\n********\n ADDED DATA TO DF RESULTS')

                    else:
                        print('DWD Station is near but not enough data')

                else:
                    print('\n********\n DWD station is not near')

            else:
                print('\n*****\n neatmo station has not data for this year')
                continue

        except Exception as msg:
            print('error while finding neighbours ', msg)
        alls_stns_len -= 1

    # drop all netatmo stns with no data
    df_results.dropna(axis=0, how='any', inplace=True)
    df_results_correlations.dropna(axis=0, how='any',
                                   inplace=True)

    # save both dataframes
    df_results.to_csv(
        os.path.join(
            out_save_dir_orig,
            'year_%s_df_comparing_p%d_and_mean_freq_%s_'
            'dwd_netatmo_with_temp_thr_%s.csv'
            % (year_val, df_min_ppt_thr,
               temp_freq_resample, use_temp_thr)),
        sep=';')
    df_results_correlations.to_csv(
        os.path.join(
            out_save_dir_orig,
            'year_%s_df_comparing_correlations_p%d_freq_%s_'
            'dwd_netatmo_with_temp_thr_%s.csv'
            % (year_val, df_min_ppt_thr,
                temp_freq_resample, use_temp_thr)),
        sep=';')

    return df_results, df_results_correlations


#==============================================================================
#
#==============================================================================


def save_how_many_abv_same_below(
    df_p1_mean,  # df with netatmo stns and result of comparing
    temp_freq,  # temp freq of df
    ppt_thr,  # used for calculating p1 or p5
    use_temp_thr,  # it True use remove all data where temp<thr
    out_dir,  # out save dir for plots
    year_vals  # if all years or year by year
):
    '''
    use this funcion to find how many stations have similar values for p1 and
    for the average, how many are netatmo are below dwd and how many 
    netatmo are above dwd

    return the result in a dataframe
    '''
    df_out = pd.DataFrame()

    df_out.loc['count_vals_same_p1',
               'Result'] = df_p1_mean[df_p1_mean['p1'] == 's'].shape[0]
    df_out.loc['count_vals_less_p1',
               'Result'] = df_p1_mean[df_p1_mean['p1'] == '_'].shape[0]
    df_out.loc['count_vals_abv_p1',
               'Result'] = df_p1_mean[df_p1_mean['p1'] == '+'].shape[0]

    df_out.loc['count_vals_same_mean',
               'Result'] = df_p1_mean[df_p1_mean['mean_ppt'] == 's'].shape[0]
    df_out.loc['count_vals_less_mean',
               'Result'] = df_p1_mean[df_p1_mean['mean_ppt'] == '_'].shape[0]
    df_out.loc['count_vals_abv_mean',
               'Result'] = df_p1_mean[df_p1_mean['mean_ppt'] == '+'].shape[0]

    df_out.to_csv(os.path.join(
        out_dir,
        'year_%_df_similarities_%s_%dmm_use_temp_thr_%s_.csv'
        % (year_vals, temp_freq, ppt_thr, use_temp_thr)))
    pass


#==============================================================================
#
#==============================================================================

def plt_on_map_comparing_p1_ppt_mean_netatmo_dwd(
    df_p1_mean,  # df with netatmo stns and result of comparing
    shp_de_file,  # shapefile of BW
    temp_freq,  # temp freq of df
    ppt_thr,  # min ppt df, select all vals abv thr
    use_temp_thr,  # it True use remove all data where temp<thr
    out_dir,  # out save dir for plots
    year_vals  # if all years or year by year
):
    '''
    Read the df_results containing for every netatmo station
    the coordinates (lon, lat) and the comparision between 
    the p1 and the mean value, between netatmo and nearest dwd

    plot the reults on a map, either with or without temp_thr
    use the shapefile of BW
    '''
    if use_temp_thr:
        title_add = 'With_Temperatre_threshold'
    else:
        title_add = '_'
    #==========================================================================
    # Plot comparing p1
    #==========================================================================
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
    for i in range(df_p1_mean.shape[0]):
        ax.scatter(df_p1_mean.lon.values[i],
                   df_p1_mean.lat.values[i],
                   marker=df_p1_mean.p1.values[i],
                   c=df_p1_mean.colorp1.values[i],
                   alpha=1,
                   s=15,
                   label='P1')

    ax.set_title('Difference in Probability P1 (Ppt>%dmm)'
                 ' Netatmo and DWD %s data year %s'
                 % (ppt_thr, temp_freq, year_vals))
    ax.grid(alpha=0.5)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect(1.0)
    plt.savefig(
        os.path.join(
            out_dir,
            'year_%s_%s_%s_differece_in_p%d_netatmo_ppt_dwd_station.png'
            % (year_vals, title_add, temp_freq, ppt_thr)),

        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)
    plt.close()
    #==========================================================================
    # Plot comparision of mean ppt
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

    for i in range(df_p1_mean.shape[0]):

        ax.scatter(df_p1_mean.lon.values[i],
                   df_p1_mean.lat.values[i],
                   marker=df_p1_mean.loc[:, 'mean_ppt'].values[i],
                   c=df_p1_mean.color_mean_ppt.values[i],
                   s=15,
                   alpha=1,
                   label='Mean')
    ax.grid(alpha=0.5)
    ax.set_title('Difference in Average Rainfall values'
                 ' Netatmo and DWD %s data year %s'
                 % (temp_freq, year_vals))
    ax.set_aspect(1.0)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    plt.savefig(
        os.path.join(
            out_dir,
            'year_%s_%s_%s_ppt_thr_%d_difference_in_mean_station_ppt_dwd_station.png'
            % (year_vals, title_add, temp_freq, ppt_thr)),
        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)
    plt.close()

    return df_p1_mean

#==============================================================================
#
#==============================================================================


def plt_on_map_agreements(
    df_correlations,  # df with netatmo stns, result of correlations
    col_to_plot,  # which column to plot , str_label_of_col
    shp_de_file,  # shapefile of BW
    temp_freq,  # temp freq of df
    ppt_thr,  # min ppt df, select all vals abv thr
    use_temp_thr,  # it True use remove all data where temp<thr
    out_dir,  # out save dir for plots
    year_vals  # if all years or year by year
):
    '''
    Read the df_results containing for every netatmo station
    the coordinates (lon, lat) and the comparision between 
    the pearson and spearman correlations,
    between netatmo and nearest dwd station
    plot the reults on a map, either with or with temp_thr
    '''
    if use_temp_thr:
        title_add = 'With_Temperatre_threshold'
    else:
        title_add = '_'

    print('plotting comparing %s' % col_to_plot)
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

    # plot the stations in shapefile, look at the results of agreements
    texts = []
    for i in range(df_correlations.shape[0]):
        ax.scatter(df_correlations.lon.values[i],
                   df_correlations.lat.values[i],
                   alpha=1,
                   c='b',
                   s=15,
                   label=df_correlations[col_to_plot].values[i])
        texts.append(
            ax.text(
                df_correlations.lon.values[i],
                df_correlations.lat.values[i],
                str(df_correlations[col_to_plot].values[i]),
                color='k'))

    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='->', color='red', lw=0.25))
    ax.set_title('%s  (Ppt>%dmm)'
                 ' Netatmo and DWD %s data year %s'
                 % (col_to_plot, ppt_thr, temp_freq, year_vals))
    ax.grid(alpha=0.5)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect(1.0)
    plt.savefig(
        os.path.join(
            out_dir,
            'year_%s_%s_%s_%s_p%d_netatmo_ppt_dwd_station.png'
            % (year_vals, title_add, temp_freq, col_to_plot, ppt_thr)),
        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)
    plt.close()
    return df_correlations

#==============================================================================
#
#==============================================================================


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    for temp_freq in aggregation_frequencies:
        print('\n********\n Time aggregation is', temp_freq)

        # call this function to get the two dfs, one containing
        # df_results of comparing p1 (or p5) and average ppt
        # df_correlations comparing correaltions, agreements
#         (df_results,
#             df_results_correlations) = compare_netatmo_dwd_p1_or_p5_or_mean_ppt(
#                 netatmo_ppt_df_file=path_to_ppt_netatmo_data,
#                 path_to_dwd_data=path_to_ppt_hdf_data,
#                 netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
#                 distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
#                 min_dist_thr_ppt=min_dist_thr_ppt,
#                 temp_freq_resample=temp_freq,
#                 df_min_ppt_thr=ppt_min_thr)
#
#         # use this funtion to find how many stations are similar,
#         # or above, or below neighbouring dwd station
#         save_how_many_abv_same_below(
#             df_p1_mean=df_results,
#             temp_freq=temp_freq,
#             ppt_thr=ppt_min_thr,
#             use_temp_thr=use_temp_thr,
#             out_dir=out_save_dir_orig)
#
#         # plot the results of df_results
#         plt_on_map_comparing_p1_ppt_mean_netatmo_dwd(
#             df_p1_mean=df_results,
#             shp_de_file=path_to_shpfile,
#             temp_freq=temp_freq,
#             ppt_thr=ppt_min_thr,
#             use_temp_thr=use_temp_thr,
#             out_dir=out_save_dir_orig, year_vals='all years')
#
#         for col_label in df_results_correlations.columns:
#             if 'Correlation' in col_label:
#                 # plot the results of df_results_correlations
#                 plt_on_map_agreements(
#                     df_correlations=df_results_correlations,
#                     col_to_plot=col_label,
#                     shp_de_file=path_to_shpfile,
#                     temp_freq=temp_freq,
#                     ppt_thr=ppt_min_thr,
#                     use_temp_thr=use_temp_thr,
#                     out_dir=out_save_dir_orig, year_vals='all years')

        for year in ['2017', '2016', '2017', '2018', '2019']:
            print('\n***\n year is', year)
            (df_results_yearly,
             df_results_correlations_yearly) = compare_netatmo_dwd_p1_or_p5_or_mean_ppt_year_by_year(
                netatmo_ppt_df_file=path_to_ppt_netatmo_data,
                path_to_dwd_data=path_to_ppt_hdf_data,
                netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
                distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
                min_dist_thr_ppt=min_dist_thr_ppt,
                temp_freq_resample=temp_freq,
                df_min_ppt_thr=ppt_min_thr,
                year_val=year)

            # plot the results of df_results
            plt_on_map_comparing_p1_ppt_mean_netatmo_dwd(
                df_p1_mean=df_results_yearly,
                shp_de_file=path_to_shpfile,
                temp_freq=temp_freq,
                ppt_thr=ppt_min_thr,
                use_temp_thr=use_temp_thr,
                out_dir=out_save_dir_orig,
                year_vals=year)

            for col_label in df_results_correlations_yearly.columns:
                if 'Correlation' in col_label:
                    # plot the results of df_results_correlations
                    plt_on_map_agreements(
                        df_correlations=df_results_correlations_yearly,
                        col_to_plot=col_label,
                        shp_de_file=path_to_shpfile,
                        temp_freq=temp_freq,
                        ppt_thr=ppt_min_thr,
                        use_temp_thr=use_temp_thr,
                        out_dir=out_save_dir_orig,
                        year_vals=year)
    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
