# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    Calculate and plot statistical differences between neighbours
Purpose: Find validity of Netatmo Station compared to DWD station

Created on: 2019-07-16

For every Netatmo precipitation station select the
convective season period (Mai till Ocotber), find the nearest
DWD station, intersect both stations, and calcualte for the 
staions the values of P1 and the mean and the Spearman and Pearson 
coorealtions and compare the two (if above, below or equal (+-10 %))
save and plot the results on a map,
see if the Netatmo stations behave similar or not.

Do it on using all data for a station

Parameters
----------

Input Files
    DWD precipitation station data
    Netatmo precipitation station data
    Netatmo station coordinates data
    Distance matrix between DWD and Netatmo stations
    Shapefile of BW area

Returns
-------
Df_results: df containing for every Netamo station, the statistical difference
    in terms of P1 (1-Prob(ppt<=1) and the average value to the nearest DWD
    station. Only stations with available neighoubrs are included
    
    
Df_table: df containing the number of Netatmo Stations with similar,
    bigger or smaller P1 and Mean values
    
Plot everything on a map using a shapefile boundaries
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


from _00_additional_functions import (calculate_probab_ppt_below_thr,
                                      resample_intersect_2_dfs,
                                      select_convective_season,
                                      compare_p1_dwd_p1_netatmo,
                                      calc_plot_contingency_table_2_stations,
                                      save_how_many_abv_same_below,
                                      plot_subplot_fig)


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
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.csv')
assert os.path.exists(path_to_ppt_netatmo_data_csv), 'wrong NETATMO Ppt file'

# for reading ppt data station by station
path_to_ppt_netatmo_data_feather = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\ppt_all_netatmo_hourly_stns_combined_.fk')
assert os.path.exists(
    path_to_ppt_netatmo_data_feather), 'wrong NETATMO Ppt file'
# for getting neighboring DWD stations
path_to_ppt_dwd_data = (
    r"F:\download_DWD_data_recent\all_dwd_hourly_ppt_data_combined_1995_2019.fk")
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'

distance_matrix_netatmo_dwd_df_file = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_DWD.csv')
assert os.path.exists(
    distance_matrix_netatmo_dwd_df_file), 'wrong Distance MTX  file'

path_to_netatmo_coords_df_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\rain_bw_1hour"
    r"\netatmo_bw_1hour_coords_with_duplicates.csv")  # TODO: CHANGE
assert os.path.exists(path_to_netatmo_coords_df_file), 'wrong DWD coords file'

path_to_netatmo_gd_stns_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\filter_Netamo_data_basedon_indicator_correlation"
    r"\keep_stns_all_neighbors_combined_95_per_.csv")
assert os.path.exists(path_to_netatmo_gd_stns_file), 'wrong netatmo stns file'

path_to_shpfile = (r'F:\data_from_exchange\Netatmo'
                   r'\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp')

assert os.path.exists(path_to_shpfile), 'wrong shapefile path'

out_save_dir_orig = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                     r'\plots_NetAtmo_ppt_DWD_ppt_neighbors_mean_and_pvals')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)


#==============================================================================
#
#==============================================================================
# used in Netatmo coords df
x_col_name = ' lon'
y_col_name = ' lat'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 500000  # 5000  # m

# threshold for max ppt value per hour
max_ppt_thr = 100.  # ppt above this value are not considered
ppt_min_thr_lst = [2]  # , 5, 10 used when calculating p1 = 1-p0

min_req_ppt_vals = 10  # minimum values that a station should have

# '120min', '480min', '720min', '1440min', 'M']
aggregation_frequencies = ['60min']

neighbors_to_chose_lst = [0, 1]  # refers to DWD neighbot (0=first)

# this is used to keep only data where month is not in this list
not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april

plot_contingency_maps_all_stns = False
plot_figures = True

date_fmt = '%Y-%m-%d %H:%M:%S'

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
        df_min_ppt_thr,  # ppt_thr, select all values above thr
        min_req_ppt_vals,  # threshold of minimum required precipitation values
        flag_plot_contingency_maps,  # if True, plot contingency maps 2 stns
):
    '''
     Find then the netatmo station the neares DWD station
     intersect both stations, calculate the difference in the
     probability that P > ppt_thr and in the mean value

     Add the result to a new dataframe and return it

    #TODO: add documenation about parameters
    '''
    print('\n######\n getting all station names, reading dfs \n#######\n')

    # get all station names for netatmo
    stns_ppt_ids = pd.read_csv(
        pth_to_netatmo_cols_df_csv, nrows=0, sep=';', engine='c',
        memory_map=True).columns.tolist()
    try:
        stns_ppt_ids = list(filter(lambda x: x != 'Unnamed: 0', stns_ppt_ids))
    except Exception as msg:
        print(msg)

    # read distance matrix dwd-netamot ppt
    in_df_distance_netatmo_dwd = pd.read_csv(
        distance_matrix_netatmo_ppt_dwd_ppt, sep=';', index_col=0)

    # read netatmo ppt coords df (for plotting)
    in_netatmo_df_coords = pd.read_csv(netatmo_ppt_coords_df, sep=',',
                                       index_col=0, engine='c')

    # read netatmo good stns df (after being filtered)
    # in_df_stns = pd.read_csv(path_netatmo_gd_stns, index_col=0,
    #                         sep=';', header=None)
    # good_stns = list(in_df_stns.values.ravel())

    print('\n######\n done reading all dfs \n#######\n')
    # create df to append results of comparing two stns
    df_results = pd.DataFrame(index=stns_ppt_ids)

    if plot_contingency_maps_all_stns:
        # append all stns with lower simultaneous values than dwd
        df_contingency_table = pd.DataFrame(index=stns_ppt_ids)

    alls_stns_len = len(stns_ppt_ids)
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

                print('\n********\n Second DWD Stn Id is', stn_2_dwd,
                      'distance is ', min_dist_ppt_dwd)

                # intersect dwd and netatmo ppt data
                if temp_freq_resample != '60min':
                    df_netatmo_cmn, df_dwd_cmn = resample_intersect_2_dfs(
                        netatmo_ppt_stn1, df_dwd, temp_freq_resample)
                else:
                    new_idx_common = netatmo_ppt_stn1.index.intersection(
                        df_dwd.index)

                    try:
                        df_netatmo_cmn = netatmo_ppt_stn1.loc[new_idx_common, :]
                        df_dwd_cmn = df_dwd.loc[new_idx_common, :]
                    except Exception:
                        df_netatmo_cmn = netatmo_ppt_stn1.loc[new_idx_common]
                        df_dwd_cmn = df_dwd.loc[new_idx_common]

                if (df_netatmo_cmn.values.shape[0] > min_req_ppt_vals and
                        df_dwd_cmn.values.shape[0] > min_req_ppt_vals):

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

                    if flag_plot_contingency_maps:
                        #======================================================
                        # construct continegency table for pairwise stations
                        #======================================================

                        df_contingency_table = calc_plot_contingency_table_2_stations(
                            first_stn_id=ppt_stn_id,  # id of first stations
                            second_stn_id=stn_2_dwd,  # id of second station
                            first_df=df_netatmo_cmn,  # df of first station
                            second_df=df_dwd_cmn,  # df of second station
                            first_thr=df_min_ppt_thr,  # threshhold for first df
                            second_thr=df_min_ppt_thr,  # threshhold for second df
                            seperating_dist=min_dist_ppt_dwd,  # distance between 2 stations
                            temp_freq=temp_freq_resample,  # temporal resolution of data
                            # df containing all stations as idx, append result
                            df_append_resutl=df_contingency_table,
                            out_plot_dir=out_save_dir_orig,  # out save dir for plot
                            plot_figures=False,  # flag is True plot all plots
                        )

                    #==========================================================
                    #  compare p1 and mean ppt for both stns
                    #==========================================================
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
                        np.nanmean(df_dwd_cmn.values),
                        np.nanmean(df_netatmo_cmn.values))

                    #==========================================================
                    # append the result to df_result, for each stn
                    #==========================================================

                    df_results.loc[ppt_stn_id,
                                   'lon'] = lon_stn_netatmo
                    df_results.loc[ppt_stn_id,
                                   'lat'] = lat_stn_netatmo
                    df_results.loc[ppt_stn_id,
                                   'Distance to neighbor'] = min_dist_ppt_dwd
                    df_results.loc[ppt_stn_id,
                                   'p%d' % df_min_ppt_thr] = compare_p1_str
                    df_results.loc[ppt_stn_id,
                                   'colorp%d' % df_min_ppt_thr] = colorp1

                    df_results.loc[ppt_stn_id, 'mean_ppt'] = compare_mean_str

                    df_results.loc[ppt_stn_id, 'color_mean_ppt'] = colormean

                    print('\n********\n ADDED DATA TO DF RESULTS')

                else:
                    print('DWD Station is near but not enough data')
            else:
                print('\n********\n DWD station is not near')

        except Exception as msg:
            print('error while finding neighbours ', msg)
            continue
    assert alls_stns_len == 0, 'not all stations were considered'

    if plot_contingency_maps_all_stns:
        df_contingency_table.dropna(axis=0, how='all', inplace=True)
        df_contingency_table.to_csv(
            os.path.join(
                out_save_dir_orig,
                'df_contingency_table'
                '_%s_ppt_thr_%0.1fmm_sep_dist_%d_neighbor_%d_.csv'
                % (temp_freq_resample,
                    df_min_ppt_thr, min_dist_thr_ppt,
                    neighbor_to_chose)))

    # save all values in one dataframes
    df_results.to_csv(
        os.path.join(out_save_dir_orig,
                     'year_allyears_df_comparing_p%d_and_mean_freq_%s_'
                     'max_sep_dist_%d_dwd_netatmo_data_considering'
                     '_neighbor_%d_.csv'
                     % (df_min_ppt_thr, temp_freq_resample, min_dist_thr_ppt,
                        neighbor_to_chose)),
        sep=';')

    return df_results


#==============================================================================
#
#==============================================================================


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    for neighbor_to_chose in neighbors_to_chose_lst:
        print('\n********\n DWD Neighbor is ', neighbor_to_chose)

        for ppt_min_thr in ppt_min_thr_lst:
            print('\n********\n Ppt Threshold is', ppt_min_thr, ' mm')

            for temp_freq in aggregation_frequencies:
                print('\n********\n Time aggregation is', temp_freq)

                # call this function to get the two dfs, one containing
                # df_results of comparing p1 (or p5) and average ppt
                # df_correlations comparing correaltions, agreements

                path_to_df_results = os.path.join(
                    out_save_dir_orig,
                    'year_allyears_df_comparing_p%d_and_mean_freq_%s_'
                    'max_sep_dist_%d_dwd_netatmo_data_considering'
                    '_neighbor_%d_.csv'
                    % (ppt_min_thr, temp_freq,
                        min_dist_thr_ppt, neighbor_to_chose))

                if (not os.path.exists(path_to_df_results) or
                        plot_contingency_maps_all_stns):

                    print('\n Data frames do not exist, creating them\n')

                    (df_results
                     ) = compare_netatmo_dwd_p1_or_p5_or_mean_ppt_or_correlations(
                        path_netatmo_ppt_df_feather=path_to_ppt_netatmo_data_feather,
                        pth_to_netatmo_cols_df_csv=path_to_ppt_netatmo_data_csv,
                        path_to_dwd_data=path_to_ppt_dwd_data,  # path_to_ppt_hdf_data,
                        path_netatmo_gd_stns=path_to_netatmo_gd_stns_file,
                        netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
                        neighbor_to_chose=neighbor_to_chose,
                        distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
                        min_dist_thr_ppt=min_dist_thr_ppt,
                        temp_freq_resample=temp_freq,
                        df_min_ppt_thr=ppt_min_thr,
                        min_req_ppt_vals=min_req_ppt_vals,
                        flag_plot_contingency_maps=plot_contingency_maps_all_stns)
                else:
                    print('\n Data frames do exist, reading them\n')
                    df_results = pd.read_csv(path_to_df_results,
                                             sep=';',
                                             index_col=0)

                if plot_figures:

                    # use these plots to find how many stations are similar,
                    # or above, or below neighbouring dwd station
                    df_table = save_how_many_abv_same_below(
                        df_p1_mean=df_results,
                        temp_freq=temp_freq,
                        ppt_thr=ppt_min_thr,
                        out_dir=out_save_dir_orig,
                        year_vals='all_years'
                    )

                    #  plot the reults on a map
                    print('\n********\n Plotting Prob maps')
                    plot_subplot_fig(
                        df_to_plot=df_results,
                        col_var_to_plot='p%d' % ppt_min_thr,
                        col_color_of_marker='colorp%d' % ppt_min_thr,
                        plot_title=('Difference in Probability P%d (Ppt>%dmm)'
                                    ' Netatmo and DWD %s data year %s'
                                    % (ppt_min_thr, ppt_min_thr,
                                       temp_freq, 'all_years')),
                        lon_col_name='lon',
                        lat_col_name='lat',
                        shapefile_path=path_to_shpfile,
                        table_to_plot=df_table,
                        out_save_name=('year_%s_%s_differece_in_p%d_netatmo'
                                       '_ppt_dwd_station.png'
                                       % ('all_years', temp_freq, ppt_min_thr)),
                        out_dir=out_save_dir_orig
                    )

                    print('\n********\n Plotting Mean maps')
                    plot_subplot_fig(
                        df_to_plot=df_results,
                        col_var_to_plot='mean_ppt',
                        col_color_of_marker='color_mean_ppt',
                        plot_title=('Difference in Average Rainfall values'
                                    ' Netatmo and DWD %s data year %s'
                                    % (temp_freq, 'all_years')),
                        lon_col_name='lon',
                        lat_col_name='lat',
                        shapefile_path=path_to_shpfile,
                        table_to_plot=df_table,
                        out_save_name=('year_%s_%s_difference'
                                       '_in_mean_station_ppt_dwd_station.png'
                                       % ('all_years', temp_freq)),
                        out_dir=out_save_dir_orig
                    )

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
