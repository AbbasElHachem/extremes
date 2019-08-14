# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Name:    ADD SCRIPT MAIN NAME
Purpose: ADD SCRIPT PURPSOSE

Created on: 2019-07-30

Parameters
----------
file_loc : str
    The file location of the spreadsheet
print_cols : bool, optional
    A flag used to print the columns to the console (default is False)

Returns
-------
list
    a list of strings representing the header columns
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

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from b_get_data import HDF5

from _00_additional_functions import (convert_coords_fr_wgs84_to_utm32_,
                                      select_df_within_period,
                                      select_convective_season,
                                      resample_intersect_2_dfs,
                                      get_cdf_part_abv_thr,
                                      plt_correlation_with_distance)
from _10_aggregate_plot_compare_2_DWD_stns import (get_dwd_stns_coords,
                                                   get_nearest_dwd_station)

from _20_claculate_statisitcal_differences_neighbouring_Netatmo_DWD_stations import (
    plt_on_map_agreements)

# =============================================================================
main_dir = Path(os.getcwd())
os.chdir(main_dir)

path_to_ppt_dwd_data = (
    r"F:\download_DWD_data_recent\all_dwd_hourly_ppt_data_combined_1995_2019.fk")
assert os.path.exists(path_to_ppt_dwd_data), 'wrong DWD Csv Ppt file'


path_to_ppt_hdf_data = (
    r'F:\download_DWD_data_recent'
    r'\DWD_60Min_ppt_stns_19950101000000_20190715000000_new.h5')

assert os.path.exists(path_to_ppt_hdf_data), 'wrong DWD Ppt file'

coords_df_file = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW"
                  r"\station_coordinates_names_hourly_only_in_BW_utm32.csv")
assert os.path.exists(coords_df_file), 'wrong DWD coords file'

path_to_shpfile = (r'F:\data_from_exchange\Netatmo'
                   r'\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp')

assert os.path.exists(path_to_shpfile), 'wrong shapefile path'


out_save_dir_orig = (
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_')


if not os.path.exists(out_save_dir_orig):
    os.mkdir(out_save_dir_orig)

#==============================================================================
#
#==============================================================================
# def epsg wgs84 and utm32 for coordinates conversion
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

x_col_name = 'X'
y_col_name = 'Y'

neighbor_to_chose = 5  # 1 refers to first neighbor
val_thr_percent = 95
aggregation_frequencies = ['60min']

not_convective_season = [10, 11, 12, 1, 2, 3, 4]


start_date = '2014-01-01 00:00:00'
end_date = '2019-07-01 00:00:00'

date_fmt = '%Y-%m-%d %H:%M:%S'

plt_figures = False  # if true plot correlations seperatly and on map
#==============================================================================
#
#==============================================================================


# @profile
def calc_indicator_correlatione_two_dwd_stns(stns_ids, tem_freq):

    in_coords_df, _, _, _ = get_dwd_stns_coords(
        coords_df_file, x_col_name, y_col_name, index_col=0,
        sep_type=';')

    stndwd_ix = ['0' * (5 - len(str(stn_id))) + str(stn_id)
                 if len(str(stn_id)) < 5 else str(stn_id)
                 for stn_id in in_coords_df.index]

    in_coords_df.index = stndwd_ix

    # intersect coordinates and stns list, get only stns in BW

    df_results_correlations = pd.DataFrame(index=stns_ids)

    stns_bw = df_results_correlations.index.intersection(in_coords_df.index)

    in_coords_df_bw = in_coords_df.loc[stns_bw, :]
    df_results_correlations_bw = df_results_correlations.loc[stns_bw, :]

    alls_stns_len = len(stns_bw)
    all_distances = []
    for iid in stns_bw:
        print('\n********\n Total number of Netatmo stations is\n********\n',
              alls_stns_len)
        alls_stns_len -= 1

        print('First Stn Id is', iid)
        try:
            #idf1 = HDF52.get_pandas_dataframe(ids=[iid])
            idf1 = pd.read_feather(path_to_ppt_dwd_data,
                                   columns=['Time', iid],
                                   use_threads=True)
            idf1.set_index('Time', inplace=True)

            idf1.index = pd.to_datetime(
                idf1.index, format=date_fmt)

            idf1.dropna(axis=0, inplace=True)

            # select only convective season
            idf1 = select_convective_season(idf1, not_convective_season)
            # select df recent years
            idf1 = select_df_within_period(idf1, start_date,
                                           end_date)
            _, stn_near, distance_near = get_nearest_dwd_station(
                first_stn_id=iid,
                coords_df_file=coords_df_file,
                x_col_name=x_col_name,
                y_col_name=y_col_name,
                index_col=0,
                sep_type=';',
                neighbor_to_chose=neighbor_to_chose)

            assert iid != stn_near, 'wrong neighbour selected'

            if len(stn_near) < 5:
                stn_near = '0' * (5 - len(str(stn_near))) + str(stn_near)
            try:
                idf2 = HDF52.get_pandas_dataframe(ids=stn_near)
            except Exception:
                raise Exception
            print('Second Stn Id is', stn_near, 'distance is', distance_near)
            all_distances.append(distance_near)

#             df_common1, df_common2 = resample_intersect_2_dfs(idf1,
#                                                               idf2,
#                                                               tem_freq)
            # intersect dwd and netatmo ppt data
            if tem_freq != '60min':
                df_common1, df_common2 = resample_intersect_2_dfs(
                    idf1, idf2, tem_freq)
            else:
                new_idx_common = idf1.index.intersection(
                    idf2.index)

                try:
                    df_common1 = idf1.loc[new_idx_common, :]
                    df_common2 = idf2.loc[new_idx_common, :]
                except Exception:
                    df_common1 = idf1.loc[new_idx_common]
                    df_common2 = idf2.loc[new_idx_common]

            if (df_common1.values.shape[0] > 0 and
                    df_common2.values.shape[0] > 0):
                df_common1 = pd.DataFrame(
                    data=df_common1.values,
                    index=df_common1.index,
                    columns=[iid])

                df_common2 = pd.DataFrame(
                    data=df_common2.values,
                    index=df_common2.index,
                    columns=[stn_near])
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
                    pears(df_common1.values,
                          df_common2.values)[0], 2)

                orig_spr_corr = np.round(
                    spr(df_common1.values,
                        df_common2.values)[0], 2)

                #==========================================================
                # select only upper tail of values of both dataframes
                #==========================================================
                val_thr_float = val_thr_percent / 100

                dwd1_cdf_x, dwd1_cdf_y = get_cdf_part_abv_thr(
                    df_common1.values, -0.1)

                # get dwd1 ppt thr from cdf
                dwd1_ppt_thr_per = dwd1_cdf_x[np.where(
                    dwd1_cdf_y >= val_thr_float)][0]

                dwd2_cdf_x, dwd2_cdf_y = get_cdf_part_abv_thr(
                    df_common2.values, -0.1)

                # get dwd2 ppt thr from cdf
                dwd2_ppt_thr_per = dwd2_cdf_x[np.where(
                    dwd2_cdf_y >= val_thr_float)][0]

                print('\n****transform values to booleans*****\n')

                df_dwd1_cmn_Bool = (
                    df_common1 > dwd1_ppt_thr_per).astype(int)
                df_dwd2_cmn_Bool = (
                    df_common2 > dwd2_ppt_thr_per).astype(int)

                # calculate spearman correlations of booleans 1, 0

                bool_spr_corr = np.round(
                    spr(df_dwd1_cmn_Bool.values.ravel(),
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
                    'Orig_Pearson_Correlation'] = orig_pears_corr

                df_results_correlations_bw.loc[
                    iid,
                    'Orig_Spearman_Correlation'] = orig_spr_corr

                df_results_correlations_bw.loc[
                    iid,
                    'Bool_Spearman_Correlation'] = bool_spr_corr

        except Exception as msg:
            print(msg)
            continue
    # assert len(all_distances) == len(stns_bw), 'smtg went wrong'
    df_results_correlations_bw.to_csv(
        os.path.join(out_save_dir_orig,
                     'year_allyears_df_dwd_correlations'
                     'freq_%s_dwd_netatmo_upper_%d_percent_data_considered'
                     '_neighbor_%d_.csv'
                     % (tem_freq,
                        val_thr_percent, neighbor_to_chose)),
        sep=';')

    return df_results_correlations_bw


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    HDF52 = HDF5(infile=path_to_ppt_hdf_data)
    dwd_ids = HDF52.get_all_ids()

    for temp_freq in aggregation_frequencies:

        df_results_correlations = calc_indicator_correlatione_two_dwd_stns(dwd_ids,
                                                                           temp_freq)

        if plt_figures:
            print('Plotting figures')
            plt_correlation_with_distance(
                df_correlations=df_results_correlations,
                dist_col_to_plot='Distance to neighbor',
                corr_col_to_plot='Bool_Spearman_Correlation',
                temp_freq=temp_freq,
                out_dir=out_save_dir_orig,
                year_vals='all_years',
                val_thr_percent=val_thr_percent,
                neighbor_nbr=neighbor_to_chose)

            plt_on_map_agreements(df_correlations=df_results_correlations,
                                  col_to_plot='Bool_Spearman_Correlation',
                                  shp_de_file=path_to_shpfile,
                                  temp_freq=temp_freq,
                                  ppt_thr=np.nan,
                                  out_dir=out_save_dir_orig,
                                  year_vals=('all_years_neighbor_%d_'
                                             % (neighbor_to_chose)),
                                  val_thr_percent=val_thr_percent)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
