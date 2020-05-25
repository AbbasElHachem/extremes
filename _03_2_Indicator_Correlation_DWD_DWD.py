# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Calculate and plot statistical differences between neighbours
Purpose: Find validity of Netatmo Station compared to DWD station

Created on: 2019-07-16

For every DWD precipitation station select the
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
    DWD precipitation station data
    DWD station coordinates data
    Shapefile of BW area

Returns
-------

    
Df_correlations: df containing for every DWD station,
    the statistical difference in terms of Pearson and Spearman 
    Correlations for original data and boolean transfomed data
    compared with the nearest DWD station.
    
Plot and save everything on a map using shapefile boundaries
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

from pathlib import Path

import os
import timeit
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears
from scipy.spatial import cKDTree


from _00_functions import (select_df_within_period,
                           select_convective_season,
                           resample_intersect_2_dfs,
                           get_cdf_part_abv_thr,
                           plt_on_map_agreements,
                           plt_correlation_with_distance)

from _01_2_read_hdf5 import HDF5
# =============================================================================

# HOURLY DATA

path_to_ppt_dwd_data_hdf5 = (
    r"X:\staff\elhachem\2020_05_20_Netatmo_CML\dwd_comb_60min_data.h5")
assert os.path.exists(path_to_ppt_dwd_data_hdf5), 'wrong DWD Csv Ppt file'

# coords of stns in utm32 in DE

# COORDINATES
path_to_dwd_coords_utm32 = (
    r"X:\staff\elhachem\2020_05_20_Netatmo_CML\DWD_60min_metadata_utm32.csv")


out_save_dir_orig = (
    r'X:\staff\elhachem\2020_05_20_Netatmo_CML\indicator_correlation')

if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

# =============================================================================
# def epsg wgs84 and utm32 for coordinates conversion and plotting coordinates
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

x_col_name = 'X'
y_col_name = 'Y'

# only highest x% of the values are selected
lower_percentile_val_lst = [99.0]  # 80, 85, 90,

# temporal aggregation of df

aggregation_frequencies = ['60min']

# month number, no need to change
# not_convective_season = [10, 11, 12, 1, 2, 3, 4]

not_convective_season = [11, 12, 13, 1, 2, 3]

# starts with one
# , 2, 3, 4, 5]  # list of which neighbors to chose
neighbors_to_chose_lst = [1]  # 4, 5, 6, 7, 8]  # 1
max_dist_thr = 100 * 1e4  # 20km
min_req_ppt_vals = 30 * 2 * 24  # stations minimum required ppt values

# select data only within this period (same as netatmo)
start_date = '2018-01-01 00:00:00'
end_date = '2019-12-31 00:00:00'

date_fmt = '%Y-%m-%d %H:%M:%S'

plt_figures = False  # if true plot correlations seperatly and on map
#==============================================================================
#
#==============================================================================


# @profile
def calc_indicator_correlatione_two_dwd_stns(
    path_to_dwd_data_hdf5,  # path_to_dwd_data_hdf5
    path_to_dwd_coords_utm32,  # path_to_dwd_coords_utm32
    tem_freq,  # temporal frequency of dataframe
    neighbor_to_chose,  # which DWD neighbor to chose (starts with 1)
    val_thr_percent,  # probabilistic percentage threshold
    min_req_ppt_vals  # stations minimum required ppt values
):

    in_dwd_df_coords_utm32 = pd.read_csv(
        path_to_dwd_coords_utm32, index_col=0,
        sep=';')

    HDF5_DWD = HDF5(infile=path_to_dwd_data_hdf5)

    all_dwd_stns_ids = HDF5_DWD.get_all_names()

    dwd_coords_xy = [(x, y) for x, y in zip(
        in_dwd_df_coords_utm32.loc[:, 'X'].values,
        in_dwd_df_coords_utm32.loc[:, 'Y'].values)]

    # create a tree from coordinates
    dwd_points_tree = cKDTree(dwd_coords_xy)
    dwd_stns_ids = in_dwd_df_coords_utm32.index

    df_results_correlations = pd.DataFrame(index=all_dwd_stns_ids)

    alls_stns_len = len(all_dwd_stns_ids)
    all_distances = []

    for iid in all_dwd_stns_ids:

        print('\n********\n Total number of DWD stations is\n********\n',
              alls_stns_len)
        alls_stns_len -= 1

        print('First Stn Id is', iid)

        try:
            # read for DWD station
            try:
                idf1 = HDF5_DWD.get_pandas_dataframe(iid)
            except Exception as msg:
                print('error reading dwd', msg)

            idf1 = select_df_within_period(idf1,
                                           start=start_date,
                                           end=end_date)

            idf1.dropna(axis=0, inplace=True)

            if idf1.size > min_req_ppt_vals:
                # find distance to all dwd stations, sort them, select minimum
                (xdwd, ydwd) = (
                    in_dwd_df_coords_utm32.loc[iid, 'X'],
                    in_dwd_df_coords_utm32.loc[iid, 'Y'])

                distances, indices = dwd_points_tree.query(
                    np.array([xdwd, ydwd]),
                    k=5)

                stn_near = dwd_stns_ids[indices[neighbor_to_chose]]

                min_dist_ppt_dwd = np.round(distances[neighbor_to_chose], 2)

                if min_dist_ppt_dwd <= max_dist_thr:
                    print('Second Stn Id is', stn_near,
                          'distance is', min_dist_ppt_dwd)
                    all_distances.append(min_dist_ppt_dwd)
                    # read for DWD station
                    try:
                        idf2 = HDF5_DWD.get_pandas_dataframe(stn_near)
                    except Exception as msg:
                        print('error reading dwd', msg)

                    idf2 = select_df_within_period(idf2,
                                                   start=start_date,
                                                   end=end_date)

                    idf2.dropna(axis=0, inplace=True)

                    if (idf1.values.ravel().shape[0] > 1 and
                            idf2.values.ravel().shape[0] > 1):
                        try:

                            df_common1, df_common2 = resample_intersect_2_dfs(
                                idf1, idf2, tem_freq)
                        except Exception as msg:
                            raise Exception
                    print('\n done resampling data')

                    if (df_common1.values.size > min_req_ppt_vals and
                            df_common2.values.size > min_req_ppt_vals):

                        # make them a dataframe
                        df_common1 = pd.DataFrame(
                            data=df_common1.values,
                            index=df_common1.index,
                            columns=[iid])

                        df_common2 = pd.DataFrame(
                            data=df_common2.values,
                            index=df_common2.index,
                            columns=[stn_near])

                        # select season
                        df_common1_summer = select_convective_season(
                            df=df_common1,
                            month_lst=not_convective_season)

                        df_common2_summer = select_convective_season(
                            df=df_common2,
                            month_lst=not_convective_season)

                        assert df_common1_summer.size == df_common2_summer.size
                        assert (df_common1_summer.isna().sum() ==
                                df_common2_summer.isna().sum(), 'NANS IN DFS')


#                         print('enough data are available for plotting')

                        # get coordinates of dwd station for plotting
#                         x_stn_dwd = in_dwd_df_coords_utm32.loc[
#                             iid, x_col_name]
#                         y_stn_dwd = in_dwd_df_coords_utm32.loc[
#                             iid, y_col_name]
#                         # convert to lon, lat (for plotting in shapefile)
#                         lon_dwd, lat_dwd = convert_coords_fr_wgs84_to_utm32_(
#                             utm32, wgs82, x_stn_dwd, y_stn_dwd)

                        # calculate pearson and spearman between original
                        # values
                        orig_pears_corr = np.round(
                            pears(df_common1_summer.values.ravel(),
                                  df_common2_summer.values.ravel())[0], 2)

                        orig_spr_corr = np.round(
                            spr(df_common1_summer.values.ravel(),
                                df_common2_summer.values.ravel())[0], 2)

                        #======================================================
                        # select only upper tail of values of both dataframes
                        #======================================================
                        val_thr_float = val_thr_percent / 100

                        dwd1_cdf_x, dwd1_cdf_y = get_cdf_part_abv_thr(
                            df_common1_summer.values.ravel(), -0.1)

                        # get dwd1 ppt thr from cdf
                        dwd1_ppt_thr_per = dwd1_cdf_x[np.where(
                            dwd1_cdf_y >= val_thr_float)][0]

                        dwd2_cdf_x, dwd2_cdf_y = get_cdf_part_abv_thr(
                            df_common2_summer.values.ravel(), -0.1)

                        # get dwd2 ppt thr from cdf
                        dwd2_ppt_thr_per = dwd2_cdf_x[np.where(
                            dwd2_cdf_y >= val_thr_float)][0]

                        print('\n****transform values to booleans*****\n')

                        df_dwd1_cmn_Bool = (
                            df_common1_summer > dwd1_ppt_thr_per).astype(int)
                        df_dwd2_cmn_Bool = (
                            df_common2_summer > dwd2_ppt_thr_per).astype(int)

                        # calculate spearman correlations of booleans 1, 0

    #                     bool_spr_corr = np.round(
    #                         spr(df_dwd1_cmn_Bool.values.ravel(),
    #                             df_dwd2_cmn_Bool.values.ravel())[0], 2)
                        bool_pears_corr = np.round(
                            pears(df_dwd1_cmn_Bool.values.ravel(),
                                  df_dwd2_cmn_Bool.values.ravel())[0], 2)

                        if bool_pears_corr < 0.2:
                            raise Exception
                            print('CHECK WASSUP')

                        #======================================================
                        # append the result to df_correlations, for each stn
                        #======================================================
#                         df_results_correlations_bw.loc[iid,
#                                                        'lon'] = lon_dwd
#                         df_results_correlations_bw.loc[iid,
#                                                        'lat'] = lat_dwd
                        df_results_correlations.loc[
                            iid,
                            'Distance to neighbor'] = min_dist_ppt_dwd

                        df_results_correlations.loc[
                            iid,
                            'DWD neighbor ID'] = stn_near
                        df_results_correlations.loc[
                            iid,
                            'DWD_orig_%s_Per_ppt_thr'
                            % val_thr_percent] = dwd1_ppt_thr_per

                        df_results_correlations.loc[
                            iid,
                            'DWD_neighb_%s_Per_ppt_thr'
                            % val_thr_percent] = dwd2_ppt_thr_per

                        df_results_correlations.loc[
                            iid,
                            'Orig_Pearson_Correlation'] = orig_pears_corr

                        df_results_correlations.loc[
                            iid,
                            'Orig_Spearman_Correlation'] = orig_spr_corr

                        df_results_correlations.loc[
                            iid,
                            'Bool_Pearson_Correlation'] = bool_pears_corr
                    else:
                        print('not enough data')
                        continue
        except Exception as msg:
            print(msg)
            continue
    print('SAVING DF')
    # assert len(all_distances) == len(stns_bw), 'smtg went wrong'

    df_results_correlations.dropna(how='all', inplace=True)
    df_results_correlations.to_csv(
        os.path.join(out_save_dir_orig,
                     'dwd_dwd_indic_corr_sep_dist_%d_'
                     'freq_%s_upper_%s_per'
                     '_neighbor_%d_.csv'  # filtered_95
                     % (max_dist_thr, temp_freq,
                        str(val_thr_percent).replace('.', '_'),
                         neighbor_to_chose)),
        sep=';')

    return df_results_correlations


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    for lower_percentile_val in lower_percentile_val_lst:
        print('\n***** lower_percentile_val: ', lower_percentile_val)

        for temp_freq in aggregation_frequencies:
            print('\n***** Temp frequency: ', temp_freq, '****')

            for neighbor_to_chose in neighbors_to_chose_lst:
                print('\n***** DWD Neighbor is: ', neighbor_to_chose, '****')

                df_results_correlations = calc_indicator_correlatione_two_dwd_stns(
                    path_to_dwd_coords_utm32=path_to_dwd_coords_utm32,
                    path_to_dwd_data_hdf5=path_to_ppt_dwd_data_hdf5,
                    tem_freq=temp_freq,
                    neighbor_to_chose=neighbor_to_chose,
                    val_thr_percent=lower_percentile_val,
                    min_req_ppt_vals=min_req_ppt_vals)

#                 if plt_figures:
#                     print('Plotting figures')
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
#                     plt_on_map_agreements(
#                         df_correlations=df_results_correlations,
#                         col_to_plot='Bool_Spearman_Correlation',
#                         shp_de_file=path_to_shpfile,
#                         temp_freq=temp_freq,
#                         out_dir=out_save_dir_orig,
#                         year_vals=('all_years_neighbor_%d_'
#                                    % (neighbor_to_chose)),
#                         val_thr_percent=lower_percentile_val)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
