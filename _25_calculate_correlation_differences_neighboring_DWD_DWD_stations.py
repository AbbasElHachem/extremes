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

# from b_get_data import HDF5

from _00_additional_functions import (convert_coords_fr_wgs84_to_utm32_,
                                      select_df_within_period,
                                      select_convective_season,
                                      resample_intersect_2_dfs,
                                      get_cdf_part_abv_thr,
                                      get_dwd_stns_coords,
                                      get_nearest_dwd_station,
                                      plt_on_map_agreements,
                                      plt_correlation_with_distance)


# =============================================================================
main_dir = Path(os.getcwd())
os.chdir(main_dir)

# path_to_ppt_dwd_data = (
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
#     r"\all_dwd_hourly_ppt_data_combined_2015_2019_.fk")
# assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'
#
# path_to_ppt_csv_data = (
#     r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
#     r"\all_dwd_hourly_ppt_data_combined_2015_2019_.csv")
#
# assert os.path.exists(path_to_ppt_csv_data), 'wrong DWD Ppt file'
#
# coords_df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
#                   r"\station_coordinates_names_hourly_only_in_BW_utm32.csv")
# assert os.path.exists(coords_df_file), 'wrong DWD coords file'

# for RH

main_dir = Path(r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz"
                )
# r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz')

# r"X:\staff\elhachem\2020_10_03_Rheinland_Pfalz"
path_to_ppt_dwd_data = (
    main_dir /
    r"ppt_dwd_2014_2019_60min_no_freezing_5deg.fk")
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'

path_to_ppt_csv_data = (
    main_dir /
    r"ppt_dwd_2014_2019_60min_no_freezing_5deg.csv")

assert os.path.exists(path_to_ppt_csv_data), 'wrong DWD Ppt file'

coords_df_file = (main_dir /
                  r"dwd_coords_in_around_RH_utm32.csv")
assert os.path.exists(coords_df_file), 'wrong DWD coords file'

# path_to_shpfile = (r"P:\2020_DFG_Netatmo\02_WPs\02_WP2\00_shapefiles"
#                    r"\BW_Landesgrenze_WGS84_UTM32N\Landesgrenze_WGS84.shp")
#
# assert os.path.exists(path_to_shpfile), 'wrong shapefile path'

# out_save_dir_orig = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#                      r'\plots_DWD_ppt_DWD_ppt_correlation_')

# out_save_dir_orig = (r'X:\staff\elhachem\Netatmo_2020\pwsflt_testdata'
#                      r'\plots_DWD_ppt_DWD_ppt_correlation_')

out_save_dir_orig = (main_dir / r'indicator_correlation')
if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

#==============================================================================
#
#==============================================================================
# def epsg wgs84 and utm32 for coordinates conversion and plotting coordinates
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

x_col_name = 'X'
y_col_name = 'Y'

# only highest x% of the values are selected
lower_percentile_val_lst = [99.0]  # 80, 85, 90,

# temporal aggregation of df
# , '120min', '480min', '720min', '1440min']
# , '120min', '480min', '720min', '1440min']
# , '120min', '480min', '720min', '1440min']
aggregation_frequencies = ['60min']

# month number, no need to change
# not_convective_season = [10, 11, 12, 1, 2, 3, 4]

not_convective_season = []

# starts with one
# , 2, 3, 4, 5]  # list of which neighbors to chose
neighbors_to_chose_lst = [1, 2]  # 4, 5, 6, 7, 8]  # 1
max_dist_thr = 100 * 1e4  # 20km
min_req_ppt_vals = 30  # stations minimum required ppt values

# select data only within this period (same as netatmo)
start_date = '2015-01-01 00:00:00'
end_date = '2019-12-30 00:00:00'

date_fmt = '%Y-%m-%d %H:%M:%S'

plt_figures = False  # if true plot correlations seperatly and on map
#==============================================================================
#
#==============================================================================


# @profile
def calc_indicator_correlatione_two_dwd_stns(
    stns_ids,  # list of all DWD stns
    tem_freq,  # temporal frequency of dataframe
    neighbor_to_chose,  # which DWD neighbor to chose (starts with 1)
    val_thr_percent,  # probabilistic percentage threshold
    min_req_ppt_vals  # stations minimum required ppt values
):

    in_coords_df, _, _, _ = get_dwd_stns_coords(
        coords_df_file, x_col_name, y_col_name, index_col_name='ID',
        sep_type=';')

#     stndwd_ix = ['0' * (5 - len(str(stn_id))) + str(stn_id)
#                  if len(str(stn_id)) < 5 else str(stn_id)
#                  for stn_id in in_coords_df.index]

#     in_coords_df.index = stndwd_ix

    # intersect coordinates and stns list, get only stns in BW

    df_results_correlations = pd.DataFrame(index=stns_ids)

    stns_bw = df_results_correlations.index.intersection(in_coords_df.index)

    in_coords_df_bw = in_coords_df.loc[stns_bw, :]

    df_results_correlations_bw = df_results_correlations.loc[stns_bw, :]

    alls_stns_len = len(stns_bw)
    all_distances = []
    for iid in stns_bw:

        #         if iid == '02088' or iid == '07331':
        #             raise Exception

        print('\n********\n Total number of DWD stations is\n********\n',
              alls_stns_len)
        alls_stns_len -= 1

        print('First Stn Id is', iid)
        try:
            # read for DWD station
            try:
                idf1 = pd.read_feather(path_to_ppt_dwd_data,
                                       columns=[
                                           'Time', iid],
                                       use_threads=True)
                idf1.set_index('Time', inplace=True)
            except Exception as msg:
                #print('error reading dwd', msg)
                idf1 = pd.read_feather(path_to_ppt_dwd_data,
                                       columns=[
                                           'index', iid],
                                       use_threads=True)
                idf1.set_index('index', inplace=True)

            idf1.index = pd.to_datetime(
                idf1.index, format=date_fmt)

            idf1.dropna(axis=0, inplace=True)

            # select only convective season
            idf1 = select_convective_season(idf1, not_convective_season)
            # select df recent years
            idf1 = select_df_within_period(idf1,
                                           start_date,
                                           end_date)

            # get id, coordinates and distances of neighbor
            _, stn_near, distance_near = get_nearest_dwd_station(
                first_stn_id=iid,
                coords_df_file=coords_df_file,
                x_col_name=x_col_name,
                y_col_name=y_col_name,
                index_col_name='ID',
                sep_type=';',
                neighbor_to_chose=neighbor_to_chose)

            assert iid != stn_near, 'wrong neighbour selected'

            if len(stn_near) < 5:
                stn_near = '0' * (5 - len(str(stn_near))) + str(stn_near)
            try:
                #                 idf2 = HDF52.get_pandas_dataframe(ids=stn_near)
                try:
                    idf2 = pd.read_feather(path_to_ppt_dwd_data,
                                           columns=[
                                               'Time', stn_near],
                                           use_threads=True)
                    idf2.set_index('Time', inplace=True)
                except Exception as msg:
                    #print('error reading dwd', msg)
                    idf2 = pd.read_feather(path_to_ppt_dwd_data,
                                           columns=[
                                               'index', stn_near],
                                           use_threads=True)
                    idf2.set_index('index', inplace=True)

                idf2.index = pd.to_datetime(
                    idf2.index, format=date_fmt)

                idf2.dropna(axis=0, inplace=True)

                idf2 = select_convective_season(idf2, not_convective_season)

                idf2 = select_df_within_period(idf2,
                                               start_date,
                                               end_date)

            except Exception:
                raise Exception

            print('Second Stn Id is', stn_near, 'distance is', distance_near)
            if distance_near < max_dist_thr:
                all_distances.append(distance_near)

                print('\n resampling data')

                if (idf1.values.ravel().shape[0] > min_req_ppt_vals and
                        idf2.values.ravel().shape[0] > min_req_ppt_vals):
                    try:

                        df_common1, df_common2 = resample_intersect_2_dfs(
                            idf1, idf2, tem_freq)
                    except Exception as msg:
                        raise Exception
                print('\n done resampling data')

                if (df_common1.values.size > min_req_ppt_vals and
                        df_common2.values.size > min_req_ppt_vals):
                    df_common1 = pd.DataFrame(
                        data=df_common1.values,
                        index=df_common1.index,
                        columns=[iid])

                    df_common2 = pd.DataFrame(
                        data=df_common2.values,
                        index=df_common2.index,
                        columns=[stn_near])
    #
    #                 plt.ioff()
    #                 fig = plt.figure(figsize=(16, 12), dpi=150)
    #                 plt.plot(df_common1.index, df_common1.values,
    #                          c='b', alpha=0.35)
    #                 plt.plot(df_common2.index, df_common2.values,
    #                          c='r', alpha=0.35)
    #
    #                 fig.savefig(os.path.join(out_save_dir_orig,
    #                                          r'dwd_%s_dwd_%s_time_%s.png'
    #                                          % (iid, stn_near, tem_freq)),
    #                             frameon=True, papertype='a4',
    #                             bbox_inches='tight', pad_inches=.2)
    #                 plt.close()
                    print('enough data are available for plotting')

                    # get coordinates of dwd station for plotting
                    x_stn_dwd = in_coords_df_bw.loc[
                        iid, x_col_name]
                    y_stn_dwd = in_coords_df_bw.loc[
                        iid, y_col_name]
                    # convert to lon, lat (for plotting in shapefile)
                    lon_dwd, lat_dwd = convert_coords_fr_wgs84_to_utm32_(
                        utm32, wgs82, x_stn_dwd, y_stn_dwd)

                    # calculate pearson and spearman between original values
                    orig_pears_corr = np.round(
                        pears(df_common1.values.ravel(),
                              df_common2.values.ravel())[0], 2)

                    orig_spr_corr = np.round(
                        spr(df_common1.values.ravel(),
                            df_common2.values.ravel())[0], 2)

                    #==========================================================
                    # select only upper tail of values of both dataframes
                    #==========================================================
                    val_thr_float = val_thr_percent / 100

                    dwd1_cdf_x, dwd1_cdf_y = get_cdf_part_abv_thr(
                        df_common1.values.ravel(), -0.1)

                    # get dwd1 ppt thr from cdf
                    dwd1_ppt_thr_per = dwd1_cdf_x[np.where(
                        dwd1_cdf_y >= val_thr_float)][0]

                    dwd2_cdf_x, dwd2_cdf_y = get_cdf_part_abv_thr(
                        df_common2.values.ravel(), -0.1)

                    # get dwd2 ppt thr from cdf
                    dwd2_ppt_thr_per = dwd2_cdf_x[np.where(
                        dwd2_cdf_y >= val_thr_float)][0]

                    print('\n****transform values to booleans*****\n')

                    df_dwd1_cmn_Bool = (
                        df_common1 > dwd1_ppt_thr_per).astype(int)
                    df_dwd2_cmn_Bool = (
                        df_common2 > dwd2_ppt_thr_per).astype(int)

                    # calculate spearman correlations of booleans 1, 0

#                     bool_spr_corr = np.round(
#                         spr(df_dwd1_cmn_Bool.values.ravel(),
#                             df_dwd2_cmn_Bool.values.ravel())[0], 2)
                    bool_spr_corr = np.round(
                        pears(df_dwd1_cmn_Bool.values.ravel(),
                              df_dwd2_cmn_Bool.values.ravel())[0], 2)
                    #==========================================================
                    # append the result to df_correlations, for each stn
                    #==========================================================
                    df_results_correlations_bw.loc[iid,
                                                   'lon'] = lon_dwd
                    df_results_correlations_bw.loc[iid,
                                                   'lat'] = lat_dwd
                    df_results_correlations_bw.loc[
                        iid,
                        'Distance to neighbor'] = distance_near

                    df_results_correlations_bw.loc[
                        iid,
                        'DWD neighbor ID'] = stn_near
                    df_results_correlations_bw.loc[
                        iid,
                        'DWD_orig_%s_Per_ppt_thr'
                        % val_thr_percent] = dwd1_ppt_thr_per

                    df_results_correlations_bw.loc[
                        iid,
                        'DWD_neighb_%s_Per_ppt_thr'
                        % val_thr_percent] = dwd2_ppt_thr_per

                    df_results_correlations_bw.loc[
                        iid,
                        'Orig_Pearson_Correlation'] = orig_pears_corr

                    df_results_correlations_bw.loc[
                        iid,
                        'Orig_Spearman_Correlation'] = orig_spr_corr

                    df_results_correlations_bw.loc[
                        iid,
                        'Bool_Spearman_Correlation'] = bool_spr_corr
                else:
                    print('not enough data')
                    continue
        except Exception as msg:
            print(msg)
            continue
    print('SAVING DF')
    # assert len(all_distances) == len(stns_bw), 'smtg went wrong'

    df_results_correlations_bw.dropna(how='all', inplace=True)
    df_results_correlations_bw.to_csv(
        os.path.join(out_save_dir_orig,
                     'pearson_year_allyears_df_dwd_correlations'
                     'freq_%s_dwd_netatmo_upper_%s_percent_data_considered'
                     '_neighbor_%d_.csv'
                     % (tem_freq,
                        str(val_thr_percent).replace('.', '_'),
                         neighbor_to_chose)),
        sep=';')

    return df_results_correlations_bw


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    dwd_ids = pd.read_csv(
        path_to_ppt_csv_data, nrows=0, sep=';', engine='c', index_col=0,
        memory_map=True).columns.tolist()
    print('Total number of DWD stations is', len(dwd_ids), '\n' * 2)

    for lower_percentile_val in lower_percentile_val_lst:
        print('\n***** lower_percentile_val: ', lower_percentile_val)

        for temp_freq in aggregation_frequencies:
            print('\n***** Temp frequency: ', temp_freq, '****')

            for neighbor_to_chose in neighbors_to_chose_lst:
                print('\n***** DWD Neighbor is: ', neighbor_to_chose, '****')

                df_results_correlations = calc_indicator_correlatione_two_dwd_stns(
                    stns_ids=dwd_ids,
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
