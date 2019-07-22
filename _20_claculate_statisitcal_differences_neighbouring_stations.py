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


import os
import timeit
import time
import shapefile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn

from adjustText import adjust_text
from pathlib import Path
from matplotlib import rc
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters


from _00_additional_functions import (calculate_probab_ppt_below_thr,
                                      resample_intersect_2_dfs,
                                      select_convective_season,
                                      constrcut_contingency_table,
                                      compare_p1_dwd_p1_netatmo,
                                      look_agreement_df1_df2,
                                      plot_subplot_fig,
                                      select_vals_abv_percentile)

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
# used in Netatmo coords df
x_col_name = ' lon'
y_col_name = ' lat'

# min distance threshold used for selecting neighbours
min_dist_thr_ppt = 5000  # m

# threshold for max ppt value per hour
max_ppt_thr = 100.  # ppt above this value are not considered
ppt_min_thr = 1  # used when calculating p1 = 1-p0
lower_percentile_val = 80  # only highest x% of the values are selected

# aggregation_frequencies = ['60min', '480min', '1440min']
aggregation_frequencies = ['60min']


# this is used to keep only data where month is not in this list
not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april

plot_contingency_maps_all_stns = True
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
        val_thr_percent,  # value in percentage, select all values above it
        flag_plot_contingency_maps,  # if True, plot contingency maps 2 stns
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
        alls_stns_len -= 1
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

                    # select only upper tail of values
                    df_netatmo_cmn = select_vals_abv_percentile(
                        df_netatmo_cmn, val_thr_percent)

                    df_dwd_cmn = select_vals_abv_percentile(
                        df_dwd_cmn, val_thr_percent)

                    if flag_plot_contingency_maps:
                        # construct continegency table
                        (df_both_below_thr,
                         df_first_abv_second_below_thr,
                         df_first_below_second_abv_thr,
                         df_both_abv_thr) = constrcut_contingency_table(
                            ppt_stn_id,
                            stn_2_dwd,
                            df_netatmo_cmn,
                            df_dwd_cmn,
                            df_min_ppt_thr,
                            df_min_ppt_thr)

                        conf_arr = np.array([[df_both_below_thr,
                                              df_first_abv_second_below_thr],
                                             [df_first_below_second_abv_thr,
                                              df_both_abv_thr]])

                        if conf_arr.shape[0] > 0:
                            print('\n++++ plotting contingency table++++\n')
                            plt.ioff()
                            fig = plt.figure(figsize=(16, 12), dpi=100)
                            ax = fig.add_subplot(111)
                            ax2 = ax.twinx()
                            ax.set_aspect(1)

                            _ = sn.heatmap(conf_arr, annot=True,
                                           vmin=0.0, cmap=plt.get_cmap('jet'),
                                           vmax=100.0, fmt='.2f')

                            ax.set_title(
                                'Contingency table \n Data above %d percent \n'
                                'Rainfall Threshold %0.1fmm '
                                '%s vs %s \n distance: %0.1f m, time freq: %s\n'
                                % (val_thr_percent, df_min_ppt_thr,
                                   ppt_stn_id, stn_2_dwd,
                                   min_dist_ppt_dwd, temp_freq_resample))

                            ax.set_yticks([0.25, 1.25])
                            ax.set_yticklabels(
                                ['Netatmo below and DWD above',
                                    'Netatmo below and DWD below'],
                                rotation=90, color='g')
                            ax.set_xticks([])

                            ax2.set_yticks([0.25, 1.25])
                            ax2.set_yticklabels(
                                ['Netatmo above and DWD below',
                                 'Netatmo above and DWD above'],
                                rotation=-90, color='m')
                            plt.savefig(
                                os.path.join(
                                    out_save_dir_orig,
                                    'contingency_table_as_a_sequence'
                                    '_%s_%s_%s_%0.1fmm_data_abv_%dpercent_.png'
                                    % (ppt_stn_id, stn_2_dwd,
                                        temp_freq_resample,
                                       df_min_ppt_thr, val_thr_percent)))

                            plt.close(fig)
                            plt.clf()

                    # get coordinates of netatmo station for plotting
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
                    #==========================================================
                    #
                    #==========================================================
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
                    df_results_correlations.loc[ppt_stn_id,
                                                'lon'] = lon_stn_netatmo
                    df_results_correlations.loc[ppt_stn_id,
                                                'lat'] = lat_stn_netatmo
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
                        'dwd_abv_thr_p%d'
                        % df_min_ppt_thr] = events_dwd_abv_thr

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
            else:
                print('\n********\n DWD station is not near')

        except Exception as msg:
            print('error while finding neighbours ', msg)

    assert alls_stns_len == 0, 'not all stations were considered'

    # drop all netatmo stns with no data
    df_results.dropna(axis=0, how='any', inplace=True)
    df_results_correlations.dropna(axis=0, how='any', inplace=True)
    df_results_nbr_of_events.dropna(axis=0, how='any', inplace=True)

    # save all dataframes
    df_results.to_csv(
        os.path.join(out_save_dir_orig,
                     'year_allyears_df_comparing_p%d_and_mean_freq_%s_'
                     'dwd_netatmo_upper_%d_percent_data_considered.csv'
                     % (df_min_ppt_thr, temp_freq_resample, val_thr_percent)),
        sep=';')

    df_results_correlations.to_csv(
        os.path.join(out_save_dir_orig,
                     'year_allyears_df_comparing_correlations_p%d_'
                     'freq_%s_dwd_netatmo_upper_%d_percent_data_considered.csv'
                     % (df_min_ppt_thr, temp_freq_resample, val_thr_percent)),
        sep=';')

    df_results_nbr_of_events.to_csv(
        os.path.join(out_save_dir_orig,
                     'year_allyears_df_comparing_nbr_of_events_p%d_'
                     'freq_%s_dwd_netatmo_upper_%d_percent_data_considered.csv'
                     % (df_min_ppt_thr, temp_freq_resample, val_thr_percent)),
        sep=';')
    return df_results, df_results_correlations, df_results_nbr_of_events

#==============================================================================
#
#==============================================================================


def save_how_many_abv_same_below(
    df_p1_mean,  # df with netatmo stns and result of comparing
    temp_freq,  # temporal freq of df
    ppt_thr,  # used for calculating p1 or p5
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

    df_out.loc['count_vals_same_p%d' % ppt_thr,
               'Result'] = df_p1_mean[df_p1_mean['p%d' % ppt_thr] == 's'].shape[0]
    df_out.loc['count_vals_less_p%d' % ppt_thr,
               'Result'] = df_p1_mean[df_p1_mean['p%d' % ppt_thr] == '_'].shape[0]
    df_out.loc['count_vals_abv_p%d' % ppt_thr,
               'Result'] = df_p1_mean[df_p1_mean['p%d' % ppt_thr] == '+'].shape[0]

    df_out.loc['count_vals_same_mean_ppt',
               'Result'] = df_p1_mean[df_p1_mean['mean_ppt'] == 's'].shape[0]
    df_out.loc['count_vals_less_mean_ppt',
               'Result'] = df_p1_mean[df_p1_mean['mean_ppt'] == '_'].shape[0]
    df_out.loc['count_vals_abv_mean_ppt',
               'Result'] = df_p1_mean[df_p1_mean['mean_ppt'] == '+'].shape[0]

#     df_out.to_csv(os.path.join(
#         out_dir,
#         'year_%s_df_similarities_%s_%dmm_.csv'
#         % (year_vals, temp_freq, ppt_thr)))
    return df_out


#==============================================================================
#
#==============================================================================


def plt_on_map_agreements(
    df_correlations,  # df with netatmo stns, result of correlations
    col_to_plot,  # which column to plot , str_label_of_col
    shp_de_file,  # shapefile of BW
    temp_freq,  # temp freq of df
    ppt_thr,  # min ppt df, select all vals abv thr
    out_dir,  # out save dir for plots
    year_vals,  # if all years or year by year
    val_thr_percent  # consider all values above it
):
    '''
    Read the df_results containing for every netatmo station
    the coordinates (lon, lat) and the comparision between 
    the pearson and spearman correlations,
    between netatmo and nearest dwd station
    plot the reults on a map, either with or with temp_thr
    '''
    if 'Bool' in col_to_plot:
        title_add = '_%dmm_' % ppt_thr
    else:
        title_add = 'all_data'

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
        texts.append(ax.text(df_correlations.lon.values[i],
                             df_correlations.lat.values[i],
                             str(df_correlations[col_to_plot].values[i]),
                             color='k'))

    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='->', color='red', lw=0.25))
    ax.set_title('%s %s Data above %dpercent'
                 ' Netatmo and DWD %s data year %s'
                 % (col_to_plot, title_add, val_thr_percent,
                     temp_freq, year_vals))
    ax.grid(alpha=0.5)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect(1.0)
    plt.savefig(
        os.path.join(
            out_dir,
            'year_%s_%s_%s_p%d_netatmo_ppt_dwd_station_above_%dpercent_.png'
            % (year_vals, temp_freq, col_to_plot, ppt_thr, val_thr_percent)),
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

        path_to_df_results = os.path.join(
            out_save_dir_orig,
            'year_allyears_df_comparing_p%d_and_mean_freq_%s_'
            'dwd_netatmo_upper_%d_percent_data_considered.csv'
            % (ppt_min_thr, temp_freq, lower_percentile_val))

        path_to_df_correlations = os.path.join(
            out_save_dir_orig,
            'year_allyears_df_comparing_correlations_p%d_'
            'freq_%s_dwd_netatmo_upper_%d_percent_data_considered.csv'
            % (ppt_min_thr, temp_freq, lower_percentile_val))

        path_to_df_nbr_of_events = os.path.join(
            out_save_dir_orig,
            'year_allyears_df_comparing_nbr_of_events_p%d_'
            'freq_%s_dwd_netatmo_upper_%d_percent_data_considered.csv'
            % (ppt_min_thr, temp_freq, lower_percentile_val))

        if (not os.path.exists(path_to_df_results) or not
                os.path.exists(path_to_df_correlations) or not
                os.path.exists(path_to_df_nbr_of_events) or
                plot_contingency_maps_all_stns):

            print('\n Data frames do not exist, creating them\n')

            (df_results,
                df_results_correlations,
                df_results_nbr_of_events) = compare_netatmo_dwd_p1_or_p5_or_mean_ppt(
                    netatmo_ppt_df_file=path_to_ppt_netatmo_data,
                    path_to_dwd_data=path_to_ppt_hdf_data,
                    netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
                    distance_matrix_netatmo_ppt_dwd_ppt=distance_matrix_netatmo_dwd_df_file,
                    min_dist_thr_ppt=min_dist_thr_ppt,
                    temp_freq_resample=temp_freq,
                    df_min_ppt_thr=ppt_min_thr,
                    val_thr_percent=lower_percentile_val,
                    flag_plot_contingency_maps=plot_contingency_maps_all_stns)
        else:
            print('\n Data frames do exist, reading them\n')
            df_results = pd.read_csv(path_to_df_results,
                                     sep=';', index_col=0)
            df_results_correlations = pd.read_csv(path_to_df_correlations,
                                                  sep=';', index_col=0)
            df_results_nbr_of_events = pd.read_csv(path_to_df_nbr_of_events,
                                                   sep=';', index_col=0)

        # use this funtion to find how many stations are similar,
        # or above, or below neighbouring dwd station
        df_table = save_how_many_abv_same_below(
            df_p1_mean=df_results,
            temp_freq=temp_freq,
            ppt_thr=ppt_min_thr,
            out_dir=out_save_dir_orig, year_vals='all_years')

        #  plot the reults on a map
        print('\n********\n Plotting Prob maps')
        plot_subplot_fig(df_to_plot=df_results,
                         col_var_to_plot='p%d' % ppt_min_thr,
                         col_color_of_marker='colorp%d' % ppt_min_thr,
                         plot_title=('Difference in Probability P%d (Ppt>%dmm)'
                                     ' Netatmo and DWD %s data year %s'
                                     % (ppt_min_thr,
                                         ppt_min_thr, temp_freq, 'all_years')),
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
        plot_subplot_fig(df_to_plot=df_results,
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
        print('\n********\n Plotting Correlation maps')
        for col_label in df_results_correlations.columns:
            if 'Correlation' in col_label:
                # plot the results of df_results_correlations
                plt_on_map_agreements(
                    df_correlations=df_results_correlations,
                    col_to_plot=col_label,
                    shp_de_file=path_to_shpfile,
                    temp_freq=temp_freq,
                    ppt_thr=ppt_min_thr,
                    out_dir=out_save_dir_orig, year_vals='all_years',
                    val_thr_percent=lower_percentile_val)

#         plt_on_map_agreements(
#             df_correlations=df_results_nbr_of_events,
#             col_to_plot='ratio_netatmo_dwd_abv_thr_p%d' % ppt_min_thr,
#             shp_de_file=path_to_shpfile,
#             temp_freq=temp_freq,
#             ppt_thr=ppt_min_thr,
#             out_dir=out_save_dir_orig, year_vals='all_years')

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
