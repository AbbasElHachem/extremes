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

Do it on a yearly basis, in this way we look for tendency in a station 

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
Df_results: df containing for every Netamo station, the statistical difference
    in terms of P1 (1-Prob(ppt<=1) and the average value to the nearest DWD
    station. Only stations with available neighoubrs are included
    
Df_correlations: df containing for every Netamo station,
    the statistical difference in terms of Pearson and Spearman 
    Correlations for original data and boolean transfomed data
    compared with the nearest DWD station.
    
Df_table: df containing the number of Netatmo Stations with similar,
    bigger or smaller P1 and Mean values
    
Plot everything on a map using shapefile boundaries
"""


__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

# TODO: CHECK IF NEEDED, IF SO THEN BREAK INTO MANY SCRIPTS
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
                                      look_agreement_df1_df2,
                                      plot_subplot_fig)

from _20_claculate_statisitcal_differences_neighbouring_stations import (
    save_how_many_abv_same_below, plt_on_map_agreements)

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


distance_matrix_df_file_ppt_temp = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_ppt_Netatmo_temp.csv')
assert os.path.exists(distance_matrix_df_file_ppt_temp), 'wrong distance df'


# path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
#                         r'\niederschlag_deutschland'
#                         r'\1993_2016_5min_merge_nan.h5')

path_to_ppt_hdf_data = (r'E:\download_DWD_data_recent'
                        r'\DWD_60Min_ppt_stns_19950101000000_20190715000000_new.h5')

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
# for NEtatmo coords df
x_col_name = 'lon'
y_col_name = 'lat'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 5000  # m

# threshold for max ppt value per hour
max_ppt_thr = 100.
ppt_min_thr = 1  # used when calculating p1 = 1-p0

aggregation_frequencies = ['60min',  '1440min']
# aggregation_frequencies = ['60min']

# this is used to keep only data where month is not in this list
not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april


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
    df_results_nbr_of_events = pd.DataFrame(index=stns_ppt_ids)

    alls_stns_len = len(stns_ppt_ids)

    for ppt_stn_id in stns_ppt_ids:

        print('\n********\n Total number of Netatmo stations is\n********\n',
              alls_stns_len)
        # iterating through netatmo ppt stations
        alls_stns_len -= 1

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
                                       'p%d' % df_min_ppt_thr] = compare_p1_str
                        df_results.loc[ppt_stn_id,
                                       'colorp%d' % df_min_ppt_thr] = colorp1

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
                            bool_pears_corr, bool_spr_corr,
                            events_dwd_abv_thr,
                            events_netatmo_abv_thr,
                            ratio_netatmo_dwd) = look_agreement_df1_df2(
                                stn_dwd_id=stn_2_dwd,
                                stn_netatmo_id=ppt_stn_id,
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
                            'Orig_Pearson_Correlation'] = orig_pears_corr

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Orig_Spearman_Correlation'] = orig_spr_coor

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Bool_Pearson_Correlation'] = bool_pears_corr

                        df_results_correlations.loc[
                            ppt_stn_id,
                            'Bool_Spearman_Correlation'] = bool_spr_corr
                    #==========================================================
                    # find events above threshold
                    #==========================================================
                        df_results_nbr_of_events.loc[ppt_stn_id,
                                                     'lon'] = lon_stn_netatmo
                        df_results_nbr_of_events.loc[ppt_stn_id,
                                                     'lat'] = lat_stn_netatmo
                        df_results_nbr_of_events.loc[
                            ppt_stn_id,
                            'dwd_abv_thr_p%d' %
                            df_min_ppt_thr] = events_dwd_abv_thr

                        df_results_nbr_of_events.loc[
                            ppt_stn_id,
                            'netatmo_abv_thr_p%d'
                            % df_min_ppt_thr] = events_netatmo_abv_thr
                        df_results_nbr_of_events.loc[
                            ppt_stn_id,
                            'ratio_netatmo_dwd_abv_thr_p%d'
                            % df_min_ppt_thr] = ratio_netatmo_dwd

                        print('\n********\n ADDED DATA TO DF RESULTS')

                    else:
                        print('DWD Station is near but not enough data')
                        continue  # test
                else:
                    print('\n********\n DWD station is not near')
                    continue  # test

            else:
                print('\n*****\n neatmo station has not data for this year')
                continue

        except Exception as msg:
            print('error while finding neighbours ', msg)

    assert alls_stns_len == 0, 'not all stations were considered'

    # drop all netatmo stns with no data
    df_results.dropna(axis=0, how='any', inplace=True)
    df_results_correlations.dropna(axis=0, how='any',
                                   inplace=True)

    # save both dataframes
    df_results.to_csv(
        os.path.join(
            out_save_dir_orig,
            'year_%s_df_comparing_p%d_and_mean_freq_%s_'
            'dwd_netatmo.csv'
            % (year_val, df_min_ppt_thr,
               temp_freq_resample)),
        sep=';')
    df_results_correlations.to_csv(
        os.path.join(
            out_save_dir_orig,
            'year_%s_df_comparing_correlations_p%d_freq_%s_'
            'dwd_netatmo.csv'
            % (year_val, df_min_ppt_thr,
                temp_freq_resample)),
        sep=';')

    df_results_nbr_of_events.to_csv(
        os.path.join(out_save_dir_orig,
                     'year_%s_df_comparing_nbr_of_events_p%d_'
                     'freq_%s_dwd_netatmo.csv'
                     % (year_val, df_min_ppt_thr, temp_freq_resample)),
        sep=';')
    return df_results, df_results_correlations, df_results_nbr_of_events


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

        for year in ['2015', '2016', '2017', '2018', '2019']:
            print('\n***\n year is', year)

            path_to_df_results_yearly = os.path.join(
                out_save_dir_orig,
                'year_%s_df_comparing_p%d_and_mean_freq_%s_dwd_netatmo.csv'
                % (year, ppt_min_thr, temp_freq))

            path_to_df_correlations_yearly = os.path.join(
                out_save_dir_orig,
                'year_%s_df_comparing_correlations_p%d_'
                'freq_%s_dwd_netatmo.csv'
                % (year, ppt_min_thr, temp_freq))

            path_to_df_nbr_of_events = os.path.join(
                out_save_dir_orig,
                'year_%s_df_comparing_nbr_of_events_p%d_'
                'freq_%s_dwd_netatmo.csv'
                % (year, ppt_min_thr, temp_freq))

            if (not os.path.exists(path_to_df_results_yearly) or not
                    os.path.exists(path_to_df_correlations_yearly) or not
                    os.path.exists(path_to_df_nbr_of_events)):

                print('\nP.S: Data frames do not exist, creating them\n')

                (df_results_yearly,
                 df_results_correlations_yearly,
                 df_results_nbr_of_events) = compare_netatmo_dwd_p1_or_p5_or_mean_ppt_year_by_year(
                    netatmo_ppt_df_file=path_to_ppt_netatmo_data,
                    path_to_dwd_data=path_to_ppt_hdf_data,
                    netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
                    distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
                    min_dist_thr_ppt=min_dist_thr_ppt,
                    temp_freq_resample=temp_freq,
                    df_min_ppt_thr=ppt_min_thr,
                    year_val=year)
            else:
                print('\nP.S: Data frames do exist, reading them\n')
                df_results_yearly = pd.read_csv(
                    path_to_df_results_yearly,
                    sep=';', index_col=0)

                df_results_correlations_yearly = pd.read_csv(
                    path_to_df_correlations_yearly,
                    sep=';', index_col=0)

            # use this funtion to find how many stations are similar,
            # or above, or below neighbouring dwd station
            df_table = save_how_many_abv_same_below(
                df_p1_mean=df_results_yearly,
                temp_freq=temp_freq,
                ppt_thr=ppt_min_thr,
                out_dir=out_save_dir_orig, year_vals=year)

    #     plot the reults on a map, either with or without temp_thr
            print('\n********\n Plotting Prob maps')
            plot_subplot_fig(df_to_plot=df_results_yearly,
                             col_var_to_plot='p1',
                             col_color_of_marker='colorp1',
                             plot_title=('Difference in Probability P1 (Ppt>%dmm)'
                                         ' Netatmo and DWD %s data year %s'
                                         % (ppt_min_thr, temp_freq, year)),
                             lon_col_name='lon',
                             lat_col_name='lat',
                             shapefile_path=path_to_shpfile,
                             table_to_plot=df_table,
                             out_save_name=('year_%s_%s_differece_in_p%d_netatmo'
                                            '_ppt_dwd_station.png'
                                            % (year, temp_freq, ppt_min_thr)),
                             out_dir=out_save_dir_orig
                             )
            print('\n********\n Plotting Mean maps')
            plot_subplot_fig(df_to_plot=df_results_yearly,
                             col_var_to_plot='mean_ppt',
                             col_color_of_marker='color_mean_ppt',
                             plot_title=('Difference in Average Rainfall values'
                                         ' Netatmo and DWD %s data year %s'
                                         % (temp_freq, year)),
                             lon_col_name='lon',
                             lat_col_name='lat',
                             shapefile_path=path_to_shpfile,
                             table_to_plot=df_table,
                             out_save_name=('year_%s_%s_difference'
                                            '_in_mean_station_ppt_dwd_station.png'
                                            % (year, temp_freq)),
                             out_dir=out_save_dir_orig
                             )

            for col_label in df_results_correlations_yearly.columns:
                if 'Correlation' in col_label:
                    # plot the results of df_results_correlations
                    plt_on_map_agreements(
                        df_correlations=df_results_correlations_yearly,
                        col_to_plot=col_label,
                        shp_de_file=path_to_shpfile,
                        temp_freq=temp_freq,
                        ppt_thr=ppt_min_thr,
                        out_dir=out_save_dir_orig,
                        year_vals=year)
#             plt_on_map_agreements(
#                 df_correlations=df_results_nbr_of_events,
#                 col_to_plot='ratio_netatmo_dwd_abv_thr_p%d' % ppt_min_thr,
#                 shp_de_file=path_to_shpfile,
#                 temp_freq=temp_freq,
#                 ppt_thr=ppt_min_thr,
#                 out_dir=out_save_dir_orig, year_vals=year)
    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
