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
    Netatmo precipitation station data
    Netatmo station coordinates data
    Distance matrix between Netatmo and Netatmo stations
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


from adjustText import adjust_text
from pathlib import Path
from matplotlib import rc
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from _00_additional_functions import (resample_intersect_2_dfs,
                                      select_convective_season,
                                      get_cdf_part_abv_thr,
                                      plt_correlation_with_distance)

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

# path_to_ppt_netatmo_data_feather = (
#     r'F:\Netatmo_5min_data'
#     r'\ppt_all_netatmo_5min_stns_combined_.fk')
# assert os.path.exists(
#     path_to_ppt_netatmo_data_feather), 'wrong NETATMO Ppt file'

distance_matrix_netatmo_netatmo_df_file = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
    r'\NetAtmo_BW\distance_mtx_in_m_NetAtmo_ppt_Netatmo_ppt.csv')
assert os.path.exists(
    distance_matrix_netatmo_netatmo_df_file), 'wrong Distance MTX  file'

path_to_netatmo_coords_df_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
    r"\rain_bw_1hour"
    r"\netatmo_bw_1hour_coords.csv")
assert os.path.exists(path_to_netatmo_coords_df_file), 'wrong DWD coords file'

path_to_netatmo_gd_stns_file = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\plots_NetAtmo_ppt_DWD_ppt_correlation_"
    r"\keep_stns_1st_neighbor_95_per_60min_.csv")
assert os.path.exists(path_to_netatmo_gd_stns_file), 'wrong netatmo stns file'

path_to_shpfile = (r'F:\data_from_exchange\Netatmo'
                   r'\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp')

assert os.path.exists(path_to_shpfile), 'wrong shapefile path'

out_save_dir_orig = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                     r'\plots_Netatmo_ppt_Netatmo_ppt_correlation_')


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


lower_percentile_val_lst = [90]  # ,95, 85, 99]  # 80, 85, 90,, 99
# only highest x% of the values are selected

# aggregation_frequencies = ['60min', '720min', '1440min']
# ['60min', '120min', '480min', '720min', '1440min']
# ,'60min', '120min', '480min', '720min', '1440min']
aggregation_frequencies = ['60min']

neighbors_to_chose_lst = [0]  # , 1, 2, 3, 4 refers to neighbor (0=first)

# this is used to keep only data where month is not in this list
not_convective_season = [10, 11, 12, 1, 2, 3, 4]  # oct till april

min_req_ppt_vals = 30  # minimum required ppt values for a station

plot_figures = True

date_fmt = '%Y-%m-%d %H:%M:%S'

#==============================================================================
#
#==============================================================================


# @profile
def compare_netatmo_dwd_p1_or_p5_or_mean_ppt_or_correlations(
        path_netatmo_ppt_df_feather,  # path to df of all netatmo ppt stations
        pth_to_netatmo_cols_df_csv,  # path to csv file, get all columns
        path_netatmo_gd_stns,  # path to netatmo stns with high rank correls
        netatmo_ppt_coords_df,  # path to netatmo ppt coords df
        neighbor_to_chose,  # which DWD station neighbor to chose
        distance_matrix_netatmo_ppt_netatmo_ppt,  # distance all netatmo-netatmo stns
        min_dist_thr_ppt,  # distance threshold when selecting dwd neigbours
        temp_freq_resample,  # temp freq to resample dfs
        val_thr_percent,  # value in percentage, select all values above it
        min_req_ppt_vals  # min required ppt values per station
):
    '''
    For every netatmo precipitation station,
     find nearest netatmo temperature station, intersect
     both stations and remove all days where temperature < 1�C

     Find then for the netatmo station the neares DWD station
     intersect both stations, calculate the difference in the
     probability that P>1mm and in the mean value

     Add the result to a new dataframe and return it

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
    in_df_distance_netatmo_netatmo = pd.read_csv(
        distance_matrix_netatmo_ppt_netatmo_ppt, sep=';', index_col=0)

    # read netatmo ppt coords df (for plotting)
    in_netatmo_df_coords = pd.read_csv(netatmo_ppt_coords_df, sep=',',
                                       index_col=0, engine='c')

    # read netatmo good stns df
    in_df_stns = pd.read_csv(path_netatmo_gd_stns, index_col=0,
                             sep=';', header=None)
    good_stns = list(in_df_stns.values.ravel())

    print('\n######\n done reading all dfs \n#######\n')
    # create df to append results of comparing two stns
    df_results_correlations = pd.DataFrame(index=stns_ppt_ids)

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

#             netatmo_ppt_stn1 = in_netatmo_ppt_stns_df.loc[:, ppt_stn_id]
            netatmo_ppt_stn1.dropna(axis=0, inplace=True)
            netatmo_ppt_stn1 = netatmo_ppt_stn1[netatmo_ppt_stn1 < max_ppt_thr]

            # select only convective season
            netatmo_ppt_stn1 = select_convective_season(
                df=netatmo_ppt_stn1,
                month_lst=not_convective_season)

            # find distance to all netatmo stations, sort them, select minimum
            distances_netatmo_to_stn1 = in_df_distance_netatmo_netatmo.loc[
                ppt_stn_id,
                in_df_distance_netatmo_netatmo.columns != ppt_stn_id]
            sorted_distances_ppt_netatmo = distances_netatmo_to_stn1.sort_values(
                ascending=True)

            # select only from neighbor to chose
            sorted_distances_ppt_netatmo = sorted_distances_ppt_netatmo.iloc[
                neighbor_to_chose:]

            # select the DWD station neighbor
            min_dist_ppt_netatmo = np.round(
                sorted_distances_ppt_netatmo.values[0], 2)

            if min_dist_ppt_netatmo <= min_dist_thr_ppt:
                # check if dwd station is near, select and read dwd stn
                stn_2_netatmo = sorted_distances_ppt_netatmo.index[0]

                df_netatmo2 = pd.read_feather(path_netatmo_ppt_df_feather,
                                              columns=['Time', stn_2_netatmo],
                                              use_threads=True)
                df_netatmo2.set_index('Time', inplace=True)

                df_netatmo2.index = pd.to_datetime(
                    df_netatmo2.index, format=date_fmt)


#                 df_dwd = HDF52.get_pandas_dataframe(ids=[stn_2_dwd])

                df_netatmo2 = df_netatmo2[df_netatmo2 < max_ppt_thr]

                # select convective season
                df_netatmo2 = select_convective_season(
                    df=df_netatmo2,
                    month_lst=not_convective_season)

                df_netatmo2.dropna(axis=0, inplace=True)
                print('\n********\n Second DWD Stn Id is', stn_2_netatmo,
                      'distance is ', min_dist_ppt_netatmo)

                # intersect dwd and netatmo ppt data
                if temp_freq_resample != '5min':
                    df_netatmo_cmn, df_netatmo2_cmn = resample_intersect_2_dfs(
                        netatmo_ppt_stn1, df_netatmo2, temp_freq_resample)
                else:
                    if (netatmo_ppt_stn1.values.shape[1] != 0 and
                            df_netatmo2.values.shape[1] != 0):

                        netatmo_ppt_stn1.dropna(how='any', inplace=True)
                        df_netatmo2.dropna(how='any', inplace=True)
                        idx_common = netatmo_ppt_stn1.index.intersection(
                            df_netatmo2.index)

                        if idx_common.shape[0] > 0:
                            print(
                                '\n********\n common index is found, interescting dataframes')
                            try:
                                df_netatmo_cmn = netatmo_ppt_stn1.loc[idx_common, :]
                                df_netatmo2_cmn = df_netatmo2.loc[idx_common, :]
                            except Exception:
                                df_netatmo_cmn = netatmo_ppt_stn1.loc[idx_common]
                                df_netatmo2_cmn = df_netatmo2.loc[idx_common]
                    else:
                        print('not enough data for intersecting 5min data')
                        continue
                if (df_netatmo_cmn.values.shape[0] > min_req_ppt_vals and
                        df_netatmo2_cmn.values.shape[0] > min_req_ppt_vals):
                    pass
                    # change everything to dataframes with stn Id as column
                    df_netatmo_cmn = pd.DataFrame(
                        data=df_netatmo_cmn.values,
                        index=df_netatmo_cmn.index,
                        columns=[ppt_stn_id])

                    df_netatmo2_cmn = pd.DataFrame(
                        data=df_netatmo2_cmn.values,
                        index=df_netatmo2_cmn.index,
                        columns=[stn_2_netatmo])

                    # get coordinates of netatmo station for plotting
                    lon_stn_netatmo = in_netatmo_df_coords.loc[
                        ppt_stn_id_name_orig, x_col_name]
                    lat_stn_netatmo = in_netatmo_df_coords.loc[
                        ppt_stn_id_name_orig, y_col_name]

                    #==========================================================
                    # look for agreements, correlation between all values
                    #==========================================================

                    # calculate pearson and spearman between original values
                    orig_pears_corr = np.round(
                        pears(df_netatmo2_cmn.values,
                              df_netatmo_cmn.values)[0], 2)

                    orig_spr_corr = np.round(
                        spr(df_netatmo2_cmn.values,
                            df_netatmo_cmn.values)[0], 2)

                    #==========================================================
                    # select only upper tail of values of both dataframes
                    #==========================================================
                    val_thr_float = val_thr_percent / 100

                    netatmo_cdf_x, netatmo_cdf_y = get_cdf_part_abv_thr(
                        df_netatmo_cmn.values.ravel(), -0.1)
                    # get netatmo ppt thr from cdf
                    netatmo_ppt_thr_per = netatmo_cdf_x[np.where(
                        netatmo_cdf_y >= val_thr_float)][0]

                    netatmo2_cdf_x, netatmo2_cdf_y = get_cdf_part_abv_thr(
                        df_netatmo2_cmn.values.ravel(), -0.1)

                    # get dwd ppt thr from cdf
                    netatmo2_ppt_thr_per = netatmo2_cdf_x[np.where(
                        netatmo2_cdf_y >= val_thr_float)][0]

                    print('\n****transform values to booleans*****\n')

                    df_netatmo_cmn_Bool = (
                        df_netatmo_cmn > netatmo_ppt_thr_per).astype(int)
                    df_netatmo2_cmn_Bool = (
                        df_netatmo2_cmn > netatmo2_ppt_thr_per).astype(int)

                    # calculate spearman correlations of booleans 1, 0

                    bool_spr_corr = np.round(
                        spr(df_netatmo2_cmn_Bool.values.ravel(),
                            df_netatmo_cmn_Bool.values.ravel())[0], 2)

                    #==========================================================
                    # append the result to df_correlations, for each stn
                    #==========================================================
                    df_results_correlations.loc[ppt_stn_id,
                                                'lon'] = lon_stn_netatmo
                    df_results_correlations.loc[ppt_stn_id,
                                                'lat'] = lat_stn_netatmo
                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Distance to neighbor'] = min_dist_ppt_netatmo

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Netatmo neighbor ID'] = stn_2_netatmo

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Netatmo_orig_%s_Per_ppt_thr'
                        % val_thr_percent] = netatmo_ppt_thr_per

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Netatmo_neigh_%s_Per_ppt_thr'
                        % val_thr_percent] = netatmo2_ppt_thr_per

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Orig_Pearson_Correlation'] = orig_pears_corr

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Orig_Spearman_Correlation'] = orig_spr_corr

                    df_results_correlations.loc[
                        ppt_stn_id,
                        'Bool_Spearman_Correlation'] = bool_spr_corr

                    print('\n********\n ADDED DATA TO DF RESULTS')

                else:
                    print('Netatmo Station is near but not enough data')
            else:
                print('\n********\n Netatmo station is not near')

        except Exception as msg:
            print('error while finding neighbours ', msg)
            continue
    assert alls_stns_len == 0, 'not all stations were considered'

    # save all dataframes
    df_results_correlations.dropna(how='all', inplace=True)
    df_results_correlations.to_csv(
        os.path.join(out_save_dir_orig,
                     'year_allyears_df_comparing_correlations_max_sep_dist_%d_'
                     'freq_%s_netatmo_netatmo_upper_%d_percent_data_considered'
                     '_neighbor_%d_.csv'
                     % (min_dist_thr_ppt, temp_freq_resample,
                        val_thr_percent, neighbor_to_chose)),
        sep=';')

    return df_results_correlations


#==============================================================================
#
#==============================================================================


def plt_on_map_agreements(
    df_correlations,  # df with netatmo stns, result of correlations
    col_to_plot,  # which column to plot , str_label_of_col
    shp_de_file,  # shapefile of BW
    temp_freq,  # temp freq of df
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
        percent_add = '_above_%dpercent' % val_thr_percent
    else:
        percent_add = '_'
    print('plotting comparing %s' % col_to_plot)
    # select values only for column and dron nans
    df_correlations = df_correlations.loc[:, [
        'lon', 'lat', col_to_plot]].dropna()

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

    ax.set_title('%s Data %s'
                 ' Netatmo and Netatmo %s data year %s'
                 % (col_to_plot,  percent_add,
                     temp_freq, year_vals))
    ax.grid(alpha=0.5)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect(1.0)
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='->', color='red', lw=0.25))

    plt.savefig(
        os.path.join(
            out_dir,
            'year_%s_%s_%s_netatmo_ppt_netatmo_station_%s.png'
            % (year_vals, temp_freq, col_to_plot, percent_add)),
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

    for lower_percentile_val in lower_percentile_val_lst:
        print('\n********\n lower_percentile_val', lower_percentile_val)

        for temp_freq in aggregation_frequencies:
            print('\n********\n Time aggregation is', temp_freq)

            for neighbor_to_chose in neighbors_to_chose_lst:
                print('\n********\n Netatmo neighbor is', neighbor_to_chose)

                # call this function to get the two dfs, one containing
                # df_results of comparing p1 (or p5) and average ppt
                # df_correlations comparing correaltions, agreements
                path_to_df_correlations = os.path.join(
                    out_save_dir_orig,
                    'year_allyears_df_comparing_correlations_max_sep_dist_%d_'
                    'freq_%s_netatmo_netatmo_upper_%d_percent_data_considered'
                    '_neighbor_%d_.csv'
                    % (min_dist_thr_ppt, temp_freq, lower_percentile_val,
                       neighbor_to_chose))

                if (not os.path.exists(path_to_df_correlations)):

                    print('\n Data frames do not exist, creating them\n')

                    (df_results_correlations
                     ) = compare_netatmo_dwd_p1_or_p5_or_mean_ppt_or_correlations(
                        path_netatmo_ppt_df_feather=path_to_ppt_netatmo_data_feather,
                        pth_to_netatmo_cols_df_csv=path_to_ppt_netatmo_data_csv,
                        path_netatmo_gd_stns=path_to_netatmo_gd_stns_file,
                        netatmo_ppt_coords_df=path_to_netatmo_coords_df_file,
                        neighbor_to_chose=neighbor_to_chose,
                        distance_matrix_netatmo_ppt_netatmo_ppt=distance_matrix_netatmo_netatmo_df_file,
                        min_dist_thr_ppt=min_dist_thr_ppt,
                        temp_freq_resample=temp_freq,
                        val_thr_percent=lower_percentile_val,
                        min_req_ppt_vals=min_req_ppt_vals)
                else:
                    print('\n Data frames do exist, reading them\n')

                    df_results_correlations = pd.read_csv(
                        path_to_df_correlations,
                        sep=';',
                        index_col=0)

                if plot_figures:

                    print('\n********\n Plotting Correlation x,y plots')
                    plt_correlation_with_distance(
                        df_correlations=df_results_correlations,
                        dist_col_to_plot='Distance to neighbor',
                        corr_col_to_plot='Bool_Spearman_Correlation',
                        temp_freq=temp_freq,
                        out_dir=out_save_dir_orig,
                        year_vals='all_years',
                        val_thr_percent=lower_percentile_val,
                        neighbor_nbr=neighbor_to_chose)

                    plt_correlation_with_distance(
                        df_correlations=df_results_correlations,
                        dist_col_to_plot='Distance to neighbor',
                        corr_col_to_plot='Orig_Spearman_Correlation',
                        temp_freq=temp_freq,
                        out_dir=out_save_dir_orig,
                        year_vals='all_years',
                        val_thr_percent=lower_percentile_val,
                        neighbor_nbr=neighbor_to_chose)

                    print('\n********\n Plotting Correlation maps')
                    for col_label in df_results_correlations.columns:
                        if ('Correlation' in col_label and
                                'Bool_Spearman' in col_label):
                            # plot the results of df_results_correlations
                            plt_on_map_agreements(
                                df_correlations=df_results_correlations,
                                col_to_plot=col_label,
                                shp_de_file=path_to_shpfile,
                                temp_freq=temp_freq,
                                out_dir=out_save_dir_orig,
                                year_vals=('all_years_%d_m_distance_neighbor_%d_'
                                           % (min_dist_thr_ppt, neighbor_to_chose)),
                                val_thr_percent=lower_percentile_val)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
