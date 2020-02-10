# !/usr/bin/env python.
# -*- coding: utf-8 -*-

#
# # from scipy.interpolate import griddata
# from osgeo import osr
#
import os
# import time
# import timeit

# import wget
# import shapefile
# import numpy as np
import wradlib as wrl
import wradlib.util as util
import osr
from matplotlib import path
# import scipy.spatial as sp
# import chardet
# import pandas as pd
# import math
# import re
# import datetime
# import os
# import pandas as pd
# import tables
# import numpy as np
# import glob
# import datetime
# import fiona
# import matplotlib.pyplot as plt
# import matplotlib
# from matplotlib.patches import Circle, Wedge, Polygon
# import os
from scipy.spatial import cKDTree
import pandas as pd
# import tables
import numpy as np
import glob
# import datetime
import fiona
import matplotlib.pyplot as plt
import matplotlib
#
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap

from adjustText import adjust_text

from _00_additional_functions import (list_all_full_path,
                                      resampleDf, select_df_within_period)


# mas for BW data
mask = np.load(r'C:\Users\hachem\Desktop\radar\masked_array_bw',
               allow_pickle=True)

#
#
main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                r'\oridinary_kriging_compare_DWD_Netatmo\Final_results')

base_path = r'C:\Users\hachem\Downloads\radolan_daily_data\*'
#base_path = r'C:\Users\hachem\Downloads\SF201907\*'
# base_path = r'C:\Users\hachem\Downloads\RW201809\*'

path_to_dwd_ppt = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                   r"\NetAtmo_BW\ppt_all_dwd_1440min_.csv")
path_to_dwd_coords = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                      r"\NetAtmo_BW\station_coordinates_names_hourly_only_in_BW.csv")
#
files = glob.glob(base_path)

ishape = fiona.open(
    r"X:\hiwi\ElHachem\Peru_Project\ancahs_dem\DEU_adm\DEU_adm1.shp")
first = ishape.next()
#
# #==============================================================================
# #
# #==============================================================================
# # DWD Coords
#
dwd_in_coords_df = pd.read_csv(path_to_dwd_coords,
                               index_col=0,
                               sep=',',
                               encoding='utf-8')
# added by Abbas, for DWD stations
stndwd_ix = ['0' * (5 - len(str(stn_id))) + str(stn_id)
             if len(str(stn_id)) < 5 else str(stn_id)
             for stn_id in dwd_in_coords_df.index]

dwd_in_coords_df.index = stndwd_ix
dwd_in_coords_df.index = list(map(str, dwd_in_coords_df.index))

dwd_lats = dwd_in_coords_df.loc[:, 'geoBreite']
dwd_lons = dwd_in_coords_df.loc[:, 'geoLaenge']

# DWD ppt
dwd_in_ppt_vals_df = pd.read_csv(
    path_to_dwd_ppt, sep=';', index_col=0, encoding='utf-8',
    engine='c')

dwd_in_ppt_vals_df.index = pd.to_datetime(
    dwd_in_ppt_vals_df.index, format='%Y-%m-%d')

dwd_in_ppt_vals_df.dropna(how='all', axis=0, inplace=True)


#==============================================================================
#
#==============================================================================
# hourly_events = [  # '2016-06-25 00:00:00',
#     '2018-06-11 14:50:00',
#     '2018-06-11 15:50:00',
#     '2018-06-11 16:50:00',
#     '2018-06-11 17:50:00',
#     '2018-06-11 18:50:00']

hourly_events = ['2018-09-06 16:50:00',
                 '2018-09-06 17:50:00',
                 '2018-09-06 18:50:00']
#'2018-09-23 17:00:00',
#'2018-09-23 18:00:00',
#'2018-09-23 19:00:00']

daily_events = [  # '2018-12-23 00:00:00',
    #'2019-05-22 00:00:00',
    '2018-05-14 00:00:00',
    '2019-07-28 00:00:00']
#
bound = [0., 1,
         2, 5, 8,
         10, 15, 20,
         25, 30, 35, 40, 45]

interval_ppt = np.linspace(0.05, 0.95)
colors_ppt = plt.get_cmap('jet_r')(interval_ppt)
cmap_ppt = LinearSegmentedColormap.from_list('name', colors_ppt)
#cmap_ppt = plt.get_cmap('jet_r')
cmap_ppt.set_over('navy')

norm = matplotlib.colors.BoundaryNorm(bound, cmap_ppt.N)


#==============================================================================
#
#==============================================================================
temp_freq = '1440min'

path_to_Qt_ok_un_first_flt__temp_flt_1st_ = main_dir / (
    r'Ppt_ok_ok_un_new2_first_flt__temp_flt__1st_%s'
    % temp_freq)
Qt_ok_un_first_flt__temp_flt_1st_ = list_all_full_path(
    '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_)

for df_file in Qt_ok_un_first_flt__temp_flt_1st_:

    if ('dwd_only') in df_file:
        print(df_file)
        path_interpolated_using_dwd = df_file

    if ('dwd_netamo') in df_file and 'unc' not in df_file:
        print(df_file)
        path_interpolated_using_netatmo_dwd = df_file

    if ('dwd_netamo') in df_file and 'unc2perc' in df_file:
        print(df_file)
        path_interpolated_using_netatmo_dwd_unc2perc = df_file

    if ('dwd_netamo') in df_file and 'unc5perc' in df_file:
        print(df_file)
        path_interpolated_using_netatmo_dwd_unc5perc = df_file

    if ('dwd_netamo') in df_file and 'unc10perc' in df_file:
        print(df_file)
        path_interpolated_using_netatmo_dwd_unc10perc = df_file

    if ('dwd_netamo') in df_file and 'unc20perc' in df_file:
        print(df_file)
        path_interpolated_using_netatmo_dwd_unc20perc = df_file

# Read interpolation
#########################################################

df_dwd = pd.read_csv(path_interpolated_using_dwd,
                     sep=';', index_col=0, parse_dates=True,
                     infer_datetime_format=True)

# f_dwd.dropna(inplace=True)

print(df_dwd.isna().sum().max())

df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                             sep=';', index_col=0,
                             parse_dates=True,
                             infer_datetime_format=True)
# df_netatmo_dwd.dropna(inplace=True)

print(df_netatmo_dwd.isna().sum().max())

df_netatmo_dwd_unc2perc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc2perc,
                                      sep=';', index_col=0,
                                      parse_dates=True, infer_datetime_format=True)

df_netatmo_dwd_unc5perc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc5perc,
                                      sep=';', index_col=0,
                                      parse_dates=True, infer_datetime_format=True)

df_netatmo_dwd_unc10perc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc10perc,
                                       sep=';', index_col=0,
                                       parse_dates=True, infer_datetime_format=True)

df_netatmo_dwd_unc20perc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc20perc,
                                       sep=';', index_col=0,
                                       parse_dates=True, infer_datetime_format=True)


#==============================================================================
#
#==============================================================================


# dwd_in_ppt_vals_df_rs = resampleDf(dwd_in_ppt_vals_df,
#                                    agg='1440min', temp_shift=-10)
#
# dwd_in_ppt_vals_df_dl = resampleDf(dwd_in_ppt_vals_df,
#                                    agg='1440min',
#                                    temp_shift=-60)

radolan_events_to_keep, event_dates = [], []
for i, file in enumerate(files):
    #     file = file + '.gz'

    event_date = ('20' + file[60:62] + '-' + file[62:64] + '-' + file[64:66]
                  + ' ' + file[66:68] + ':' + file[68:70] + ':00')

    event_date = pd.DatetimeIndex([event_date])

    event_end_date = event_date + pd.Timedelta(minutes=10)

    if event_end_date[0] in df_dwd.index:
        print('Event: ', i, '/', len(files),  event_date)
        radolan_events_to_keep.append(file)
        event_dates.append(event_end_date)


for intense_evt, event_date in zip(radolan_events_to_keep, event_dates):
    print('going through intense events')
    #     event_date = ('20' + intense_evt[60:62] + '-' + intense_evt[62:64] + '-' + intense_evt[64:66]
    #                   + ' ' + intense_evt[66:68] + ':' + intense_evt[68:70] + ':00')

    #     event_start_date = pd.datetime(year=event_date.year[0],
    #                                    month=event_date.month[0],
    #                                    day=event_date.day[0])
    #                                     hour=event_date.hour[0],
    #                                    minute=event_date.minute[0])

    #     event_start_date_dayb4 = event_start_date - pd.Timedelta(days=1)
    # accumulated sum day before
    #     df_hourly_agg_b4 = select_df_within_period(dwd_in_ppt_vals_df,
    #                                                start=event_start_date_dayb4,
    #                                                end=event_start_date).sum()
    original_ppt = dwd_in_ppt_vals_df.loc[event_date, :]
    interpolated_ppt_dwd = df_dwd.loc[event_date, :]
    interpolated_ppt_netatmo_dwd = df_netatmo_dwd.loc[event_date, :]
    interpolated_ppt_netatmo_dwd_unc2perc = df_netatmo_dwd_unc2perc.loc[
        event_date, :]
    interpolated_ppt_netatmo_dwd_unc5perc = df_netatmo_dwd_unc5perc.loc[
        event_date, :]
    interpolated_ppt_netatmo_dwd_unc10perc = df_netatmo_dwd_unc10perc.loc[
        event_date, :]
    interpolated_ppt_netatmo_dwd_unc20perc = df_netatmo_dwd_unc20perc.loc[
        event_date, :]

#         event_start_date = event_date - pd.Timedelta(hours=24, minutes=50)
#
# #         event_end_date = event_date + pd.Timedelta(minutes=10)
# #
# #         df_hourly_agg = select_df_within_period(dwd_in_ppt_vals_df,
# #                                                 start=event_start_date,
# #                                                 end=event_end_date[0]).sum()
# #         df_hourly_agg = df_hourly_agg + df_hourly_agg_b4
#         #         if event_date in hourly_events:
#         event_date_minus_one_hr = event_date - pd.Timedelta(minutes=50)
#
#         df_hourly_agg_2 = select_df_within_period(
#             dwd_in_ppt_vals_df,
#             start=event_start_date[0],
#             end=event_date_minus_one_hr[0]).sum()
#         df_hourly_agg_2 = df_hourly_agg_2 + df_hourly_agg_b4
    # df_hourly_agg_2.sum()
    # df_hourly_agg.sum()
#         event_time_rs = event_date
#         event_time = event_date + pd.Timedelta(minutes=10)
#         shifted_event = event_date_minus_one_hr + pd.Timedelta(minutes=10)

#         shifted_event_str = str(shifted_event[0]).replace('23', '00')
#         dwd_in_ppt_vals_df.loc[shifted_event_str,:]
    #     shifted_event = event_time
    rwdata, rwattrs = wrl.io.read_radolan_composite(intense_evt)
    # mask data
    sec = rwattrs['secondary']
    rwdata.flat[sec] = -9999
    rwdata = np.ma.masked_equal(rwdata, -9999)

    # create radolan projection object
    proj_stereo = wrl.georef.create_osr("dwd-radolan")

    # create wgs84 projection object
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    # get radolan grid
    radolan_grid_xy = wrl.georef.get_radolan_grid(900, 900)
    x1 = radolan_grid_xy[:, :, 0]
    y1 = radolan_grid_xy[:, :, 1]

#         # convert to lonlat
    radolan_grid_ll = wrl.georef.reproject(radolan_grid_xy,
                                           projection_source=proj_stereo,
                                           projection_target=proj_wgs)
    lon1 = radolan_grid_ll[:, :, 0]
    lat1 = radolan_grid_ll[:, :, 1]

    radolan_coords = np.array([(lo, la) for lo, la in zip(
        lon1.flatten(), lat1.flatten())])

    radolan_coords_tree = cKDTree(radolan_coords)

    dwd_in_coords_df = dwd_in_coords_df.loc[
        dwd_in_coords_df.index.intersection(original_ppt.columns), :]

    df_ppt_radolan = pd.DataFrame(index=original_ppt.columns)

    ix_lon_loc = []
    ix_lat_loc = []
    texts = []

    plt.ioff()
    fig = plt.figure(figsize=(20, 12), dpi=200)
    ax = fig.add_subplot(111)

    print('\nGetting station data\n')

    for stn, x0, y0 in zip(dwd_in_coords_df.index,
                           dwd_lons.values,
                           dwd_lats.values):

        dd, ii = radolan_coords_tree.query([x0, y0], k=6)

        grd_lon_loc = lon1.flatten()[ii]
        grd_lat_loc = lat1.flatten()[ii]
        ppt_loc = rwdata.data.flatten()[ii]

        ix_lon_loc.append(grd_lon_loc)
        ix_lat_loc.append(grd_lat_loc)
#         grd_lon_lat_coords = np.array([(lo, la) for lo, la in zip(
#             grd_lon_loc, grd_lat_loc)])

#         df_stn = dwd_in_ppt_vals_df.loc[event_time, stn]
#         df_stn2 = dwd_in_ppt_vals_df.loc[shifted_event, stn]
#         texts.append(ax.text(x0,
#                              y0,
#                              stn,
#                              color='m'))

#         plt.ioff()
#         plt.scatter(x0, y0)
#         plt.scatter(grd_lon_loc, grd_lat_loc)
#         plt.show()
#
        for i, ppt in enumerate(ppt_loc):

            df_ppt_radolan.loc[stn, 'loc_%d_1' % i] = ppt

#         plt.ioff()
#         plt.scatter(x0, y0)
#         plt.scatter(ix_lon_loc, ix_lat_loc)
#         plt.show()
    print('Done getting station data\n')
    df_ppt_radolan[df_ppt_radolan < 0] = np.nan

    df_ppt_radolan_mean_stn = df_ppt_radolan.mean(axis=1)

    # obsv
#         df_ppt_obsv_orig_rs = dwd_in_ppt_vals_df_rs.loc[
#             event_time_rs, :].dropna(how='all')
#
#         df_ppt_obsv_orig = dwd_in_ppt_vals_df_dl.loc[
#             shifted_event, :].dropna(how='all')
    # df_ppt_obsv_orig.T.plot(legend=False)
#         df_ppt_obsv_shifted = dwd_in_ppt_vals_df.loc[
#             shifted_event, :].dropna(how='all')
    # df_ppt_obsv_shifted.T.plot(legend=False)

    #======================================================================
    # # RADAR BILD
    #======================================================================
    print('Plotting radar data\n')
    plt.ioff()
#     mask = np.ones_like(lon1, dtype=np.bool)

#     for n, i_poly in enumerate(first['geometry']['coordinates']):
#
#         p = path.Path(np.array(i_poly)[0, :, :])
#         grid_mask = p.contains_points(
#             np.vstack((lon1.flatten(),
#                        lat1.flatten())).T).reshape(900, 900)
#         mask[grid_mask] = 0

    rwdata[mask] = -1
    # mask.dump(r'C:\Users\hachem\Desktop\radar\masked_array')

    rw_maskes = np.ma.masked_array(rwdata, rwdata < 0.)

#     min_x = xss[np.argmin([xs - x0 for xs in xss])]
#     min_y = yss[np.argmin([ys - y0 for ys in yss])]

    plt.pcolormesh(lon1, lat1, rw_maskes, cmap=cmap_ppt,
                   vmin=0, norm=norm, vmax=bound[-1])
    cbar = plt.colorbar(extend='max', label='[mm/h]')

    plt.scatter(ix_lon_loc, ix_lat_loc,
                marker='.', color='k', alpha=0.5, s=10)

    plt.scatter(dwd_lons.values, dwd_lats.values,
                marker='x', alpha=0.5, color='m', s=10)

#     adjust_text(texts, ax=ax,
#                 arrowprops=dict(arrowstyle='->',
#                                 color='g', lw=0.25))
    plt.title(str(event_date[0]))
    plt.ylabel('Latitude [°]')
    plt.xlabel('Longitude [°]')

    plt.xlim([7.1, 10.7])
    plt.ylim([47.3, 50.0])

#     plt.axis('equal')
    # radar.set_aspect('equal', 'box')
    # plt.xlim([netatmo.get_xlim()[0], netatmo.get_xlim()[1]])
    # plt.ylim([netatmo.get_ylim()[0], netatmo.get_ylim()[1]])

    # plt.tight_layout()
    plt.savefig(
        os.path.join(
            r'C:\Users\hachem\Desktop\radar\ppt_radar_vs_ppt_dwd',
            'radar_event_{}.png'.format(
                str(event_date[0]).split(' ')[0])))
#     plt.show()

    plt.close()
    print('Plotting station data\n')
    #======================================================================
    #
    #======================================================================

    mse_dwd_interp = np.square(
        np.nanmean(np.subtract(original_ppt.values[0],
                               interpolated_ppt_dwd.values[0])))
    mse_dwd_netatmo_interp = np.square(
        np.nanmean(np.subtract(original_ppt.values[0],
                               interpolated_ppt_netatmo_dwd.values[0])))
    mse_dwd_radar = np.square(
        np.nanmean(np.subtract(original_ppt.values[0],
                               np.nanmean(df_ppt_radolan.values))))

    interpolated_ppt_dwd.columns = original_ppt.columns
    interpolated_ppt_netatmo_dwd.columns = original_ppt.columns

    plt.figure(figsize=(20, 12), dpi=200)

    plt.plot(original_ppt.columns.to_list(),
             original_ppt.values[0], c='r', alpha=0.5, marker='+',
             label='DWD Obs')  # Obs Same Time (10min difference)

    plt.plot(interpolated_ppt_dwd.columns.to_list(),
             interpolated_ppt_dwd.values[0], c='b', alpha=0.5, marker='*',
             label='DWD Interp RMSE %0.4f' % mse_dwd_interp)

    plt.plot(interpolated_ppt_netatmo_dwd.columns.to_list(),
             interpolated_ppt_netatmo_dwd.values[0], c='m',
             alpha=0.5, marker='2',
             label='DWD-Netatmo Interp RMSE %0.4f' % mse_dwd_netatmo_interp)

    # plt.plot(df_ppt_obsv_orig.columns.to_list(),
    #        df_ppt_obsv_orig.values[0], c='b', alpha=0.5, marker='*',
    #       label='DWD stns')  # Obs Same Time (10min difference)

#         plt.plot(df_hourly_agg.index.to_list(),
#                  df_hourly_agg.values, c='r', alpha=0.5, marker='1',
#                  label='DWD stns')

#         plt.plot(df_hourly_agg_2.index.to_list(),
#                  df_hourly_agg_2.values, c='b', alpha=0.5, marker='*',
#                  label='DWD stns')
#         plt.plot(df_ppt_obsv_shifted.columns.to_list(),
#                  df_ppt_obsv_shifted.values[0], c='b', alpha=0.75,
#                  label='Obs Shifted Time (-1hr)')

    for col in df_ppt_radolan.columns:
        plt.plot(df_ppt_radolan.index,
                 df_ppt_radolan.loc[:, col].values,
                 c='k', alpha=0.25)

    plt.plot(df_ppt_radolan.index,
             df_ppt_radolan.loc[:, col].values,
             c='k', alpha=0.05, label='Radar RMSE %0.4f' % mse_dwd_radar)
#         plt.title('Radolan: %s \nDWD: %s to %s'
#                   % (str(event_date[0]),
#                      str(event_start_date[0]),
#                      str(event_date_minus_one_hr[0])))
    plt.xticks(rotation=90)
    plt.ylabel('[mm/day]')
#     plt.xlabel([df_ppt_obsv.columns.to_list()[::10]], rotation=90)
    plt.legend(loc=0)
    # plt.tight_layout()
    plt.savefig(
        os.path.join(r'C:\Users\hachem\Desktop\radar\ppt_radar_vs_ppt_dwd',
                     'ppt_event_{}.png'.format(
                         str(event_date[0]).split(' ')[0])))

    plt.close()

    print('Done for this event')
#     break
