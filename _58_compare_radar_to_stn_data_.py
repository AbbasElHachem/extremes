# !/usr/bin/env python.
# -*- coding: utf-8 -*-

import os
import osr
import glob
import fiona
import matplotlib

import wradlib as wrl
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from pathlib import Path
from scipy.spatial import cKDTree
from matplotlib.colors import LinearSegmentedColormap
# from adjustText import adjust_text

from _00_additional_functions import (list_all_full_path,
                                      get_radar_intense_events,
                                      # resampleDf,
                                      # select_df_within_period
                                      )

plt.ioff()
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'axes.labelsize': 16})
#==============================================================================
#
#==============================================================================
# mas for BW data
mask = np.load(r'C:\Users\hachem\Desktop\radar\masked_array_bw',
               allow_pickle=True)

# temp_freq = '1440min'

temp_freq = '1440min'
#
main_dir = Path(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                r'\oridinary_kriging_compare_DWD_Netatmo\Final_results')

base_path = (r'X:\exchange\ElHachem\radolan_%s_data\*' % temp_freq)
#base_path = r'C:\Users\hachem\Downloads\SF201907\*'
# base_path = r'C:\Users\hachem\Downloads\RW201809\*'

path_to_dwd_ppt = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                   r"\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq)

path_to_dwd_coords = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\NetAtmo_BW\station_coordinates_names_hourly_only_in_BW.csv")
#
files = glob.glob(base_path)

ishape = fiona.open(
    r"X:\hiwi\ElHachem\Peru_Project\ancahs_dem\DEU_adm\DEU_adm1.shp")
first = ishape.next()


out_dir = r'C:\Users\hachem\Desktop\radar\ppt_radar_vs_ppt_dwd\%s_freq' % temp_freq
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
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

#
# daily_events = [  # '2018-12-23 00:00:00',
#     #'2019-05-22 00:00:00',
#     '2018-05-14 00:00:00',
#     '2019-07-28 00:00:00']


if temp_freq == '60min':
    bound = [0., 1, 2, 5, 8, 10, 15, 20, 25, 30]
if temp_freq == '1440min':
    bound = [0., 1, 2, 5, 8, 10, 15, 20, 25, 30, 35, 40, 45]

interval_ppt = np.linspace(0.05, 0.95)
colors_ppt = plt.get_cmap('jet_r')(interval_ppt)
cmap_ppt = LinearSegmentedColormap.from_list('name', colors_ppt)
#cmap_ppt = plt.get_cmap('jet_r')
cmap_ppt.set_over('navy')

norm = matplotlib.colors.BoundaryNorm(bound, cmap_ppt.N)

#==============================================================================
#
#==============================================================================


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
# temp_freq = '1440min'
intense_events_df_index_lst = pd.read_csv(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\dwd_%s_maximum_100_event.csv"
    % temp_freq,
    index_col=0, sep=';',
    parse_dates=True,
    infer_datetime_format=True).index.to_list()
# f_dwd.dropna(inplace=True)

# print(df_dwd.isna().sum().max())

df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                             sep=';', index_col=0,
                             parse_dates=True,
                             infer_datetime_format=True)
# df_netatmo_dwd.dropna(inplace=True)

# print(df_netatmo_dwd.isna().sum().max())

df_netatmo_dwd_unc2perc = pd.read_csv(
    path_interpolated_using_netatmo_dwd_unc2perc,
    sep=';', index_col=0,
    parse_dates=True, infer_datetime_format=True)

df_netatmo_dwd_unc5perc = pd.read_csv(
    path_interpolated_using_netatmo_dwd_unc5perc, sep=';', index_col=0,
    parse_dates=True, infer_datetime_format=True)

df_netatmo_dwd_unc10perc = pd.read_csv(
    path_interpolated_using_netatmo_dwd_unc10perc, sep=';', index_col=0,
    parse_dates=True, infer_datetime_format=True)

df_netatmo_dwd_unc20perc = pd.read_csv(
    path_interpolated_using_netatmo_dwd_unc20perc, sep=';', index_col=0,
    parse_dates=True, infer_datetime_format=True)

#==============================================================================
# # radar events to keep
#==============================================================================

try:
    df_radar_events_to_keep = pd.read_csv(r'X:\exchange\ElHachem'
                                          r'\%s_intense_events_radolan_filess.csv'
                                          % temp_freq, index_col=0,
                                          sep=';')
except Exception as msg:
    df_radar_events_to_keep = get_radar_intense_events(files,
                                                       intense_events_df_index_lst,
                                                       temp_freq)


#==============================================================================
#
#==============================================================================
# dwd_in_ppt_vals_df_rs = resampleDf(dwd_in_ppt_vals_df,
#                                    agg='1440min', temp_shift=-10)
#
# dwd_in_ppt_vals_df_dl = resampleDf(dwd_in_ppt_vals_df,
#                                    agg='1440min',
#                                    temp_shift=-60)

(dwd_dwd_rmse, dwd_netatmo_rmse,
    dwd_netatmo_unc_rmse, dwd_radolan_rmse) = [], [], [], []

dwd_dwd_rmse_stns = {stn: [] for stn in df_netatmo_dwd.columns}
dwd_netatmo_rmse_stns = {stn: [] for stn in df_netatmo_dwd.columns}
dwd_netatmo_unc_rmse_stns = {stn: [] for stn in df_netatmo_dwd.columns}
dwd_radolan_rmse_stns = {stn: [] for stn in df_netatmo_dwd.columns}


# arr_tuples = [(ix, df_dwd.columns.to_list()) for ix in df_dwd.index]
# time_stns_arr = [ix, df_dwd.columns.to_list())
#                  for ix in df_dwd.index.to_list()]
# index = pd.MultiIndex.from_tuples(time_stns_arr, names=['Time', 'Stations'])

for event_date, radar_evt in zip(df_radar_events_to_keep.index,
                                 df_radar_events_to_keep.values):

    print('going through intense event: ', event_date)

    event_date = pd.DatetimeIndex([event_date])

    radolan_event_date = event_date - pd.Timedelta(minutes=10)

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

    interpolated_ppt_dwd = interpolated_ppt_dwd.reindex(
        sorted(interpolated_ppt_dwd.columns), axis=1)

    interpolated_ppt_netatmo_dwd = interpolated_ppt_netatmo_dwd.reindex(
        sorted(interpolated_ppt_netatmo_dwd.columns), axis=1)

    interpolated_ppt_netatmo_dwd_unc10perc = (
        interpolated_ppt_netatmo_dwd_unc10perc.reindex(
            sorted(interpolated_ppt_netatmo_dwd_unc10perc.columns), axis=1))

    # for plotting
    if temp_freq == '60min':
        start_minutes = 50
        event_start_date = radolan_event_date - \
            pd.Timedelta(minutes=start_minutes)
        event_date_minus_one_hr = event_date

    if temp_freq == '1440min':
        start_minutes = 50 + float(temp_freq.replace('min', ''))
        event_start_date = radolan_event_date - \
            pd.Timedelta(minutes=start_minutes)
        event_date_minus_one_hr = radolan_event_date - pd.Timedelta(minutes=50)

#     event_end_date = event_date + pd.Timedelta(minutes=10)
#
#     df_hourly_agg = select_df_within_period(dwd_in_ppt_vals_df,
#                                             start=event_start_date,
#                                             end=event_end_date[0]).sum()

    # read radolan file
    rwdata, rwattrs = wrl.io.read_radolan_composite(radar_evt[0])
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

    # convert to lonlat
    radolan_grid_ll = wrl.georef.reproject(radolan_grid_xy,
                                           projection_source=proj_stereo,
                                           projection_target=proj_wgs)
    lon1 = radolan_grid_ll[:, :, 0]
    lat1 = radolan_grid_ll[:, :, 1]
    # using the coordinates make lon, lat KDE tree, to find neighbors of stns
    radolan_coords = np.array([(lo, la) for lo, la in zip(
        lon1.flatten(), lat1.flatten())])

    radolan_coords_tree = cKDTree(radolan_coords)

    dwd_in_coords_df = dwd_in_coords_df.loc[
        dwd_in_coords_df.index.intersection(original_ppt.columns), :]

    df_ppt_radolan = pd.DataFrame(index=original_ppt.columns)

    plt.ioff()
    fig = plt.figure(figsize=(20, 20), dpi=200)
    ax = fig.add_subplot(111)

    print('\nGetting station data\n')
    ix_lon_loc, ix_lat_loc, texts = [], [], []

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

#         texts.append(ax.text(x0,
#                              y0,
#                              stn,
#                              color='m'))

        for i, ppt in enumerate(ppt_loc):
            df_ppt_radolan.loc[stn, 'loc_%d_1' % i] = ppt

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

    rwdata[mask] = -1

    rw_maskes = np.ma.masked_array(rwdata, rwdata < 0.)

    im = ax.pcolormesh(lon1, lat1, rw_maskes, cmap=cmap_ppt,
                       vmin=0, norm=norm, vmax=bound[-1])
    fig.colorbar(im, ax=ax, extend='max', label='[mm/%s]' % temp_freq)

    ax.scatter(ix_lon_loc, ix_lat_loc,
               marker='.', color='k', alpha=0.5, s=10)

    ax.scatter(dwd_lons.values, dwd_lats.values,
               marker='x', alpha=0.5, color='m', s=10)

#     adjust_text(texts, ax=ax,
#                 arrowprops=dict(arrowstyle='->',
#                                 color='g', lw=0.25))
    ax.set_title('Event date %s' % str(event_date[0]))
    ax.set_ylabel('Latitude [°]')
    ax.set_xlabel('Longitude [°]')

    ax.set_xlim([7.1, 10.7])
    ax.set_ylim([47.3, 50.0])
    ax.grid(alpha=0.5)
#     ax.axis('equal')
#     plt.tight_layout()
    plt.savefig(
        os.path.join(
            out_dir,
            'radar_event_{}.png'.format(
                str(event_date[0]).replace(' ', '_').replace(':', ''))),
        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)
#     plt.show()
    plt.close()

    #======================================================================
    #
    #======================================================================
    print('Plotting station data\n')
    # rearange columns to have same order for plotting
#     interpolated_ppt_dwd.columns = original_ppt.columns
#     interpolated_ppt_netatmo_dwd.columns = original_ppt.columns
#     interpolated_ppt_netatmo_dwd_unc10perc.columns = original_ppt.columns
    # leave common stns per all data
    cmn_stns = [stn for stn in original_ppt.columns
                if original_ppt[stn].values >= 0
                if interpolated_ppt_dwd[stn].values >= 0
                if interpolated_ppt_netatmo_dwd[stn].values >= 0
                if interpolated_ppt_netatmo_dwd_unc10perc[stn].values >= 0
                if df_ppt_radolan_mean_stn[stn] >= 0]

    original_ppt = original_ppt.loc[:, cmn_stns]
    interpolated_ppt_dwd = interpolated_ppt_dwd.loc[:, cmn_stns]
    interpolated_ppt_netatmo_dwd = interpolated_ppt_netatmo_dwd.loc[:, cmn_stns]
    interpolated_ppt_netatmo_dwd_unc10perc = (
        interpolated_ppt_netatmo_dwd_unc10perc.loc[:, cmn_stns])

    df_ppt_radolan = df_ppt_radolan.loc[cmn_stns, :]

    # calculate RMSE observed-interpoalted for this event
    mse_dwd_interp = np.nanmean(np.sqrt(np.square(
        np.subtract(original_ppt.values[0],
                    interpolated_ppt_dwd.values[0]))))
    dwd_dwd_rmse.append(mse_dwd_interp)

    mse_dwd_netatmo_interp = np.nanmean(np.sqrt(np.square(
        np.subtract(original_ppt.values[0],
                    interpolated_ppt_netatmo_dwd.values[0]))))
    dwd_netatmo_rmse.append(mse_dwd_netatmo_interp)

    mse_dwd_netatmo_unc_interp = np.nanmean(np.sqrt(np.square(
        np.subtract(original_ppt.values[0],
                    interpolated_ppt_netatmo_dwd_unc10perc.values[0]))))
    dwd_netatmo_unc_rmse.append(mse_dwd_netatmo_unc_interp)

    try:
        mse_dwd_radar = np.nanmean(np.sqrt(np.square(
            np.subtract(original_ppt.values[0],
                        np.nanmean(df_ppt_radolan_mean_stn.values)))))
    except Exception as msg:
        print(msg)
        pass
    dwd_radolan_rmse.append(mse_dwd_radar)

    # get RMSE values per stations
    for stn in cmn_stns:
        original_ppt_stn = original_ppt.loc[
            :, stn].values[0]
        interpolated_ppt_dwd_stn = interpolated_ppt_dwd.loc[
            :, stn].values[0]
        interpolated_ppt_netatmo_dwd_stn = interpolated_ppt_netatmo_dwd.loc[
            :, stn].values[0]
        interpolated_ppt_netatmo_dwd_unc_stn = (
            interpolated_ppt_netatmo_dwd_unc10perc.loc[:, stn].values[0])
        df_ppt_radolan_stn = np.nanmean(df_ppt_radolan.loc[stn, :].values)

        dwd_dwd_rmse_stns[stn].append(
            np.subtract(original_ppt_stn, interpolated_ppt_dwd_stn))
        dwd_netatmo_rmse_stns[stn].append(
            np.subtract(original_ppt_stn, interpolated_ppt_netatmo_dwd_stn))
        dwd_netatmo_unc_rmse_stns[stn].append(
            np.subtract(original_ppt_stn, interpolated_ppt_netatmo_dwd_unc_stn))
        dwd_radolan_rmse_stns[stn].append(
            np.subtract(original_ppt_stn, df_ppt_radolan_stn))

    # plot station and interpolations
    plt.figure(figsize=(24, 12), dpi=200)

    plt.plot(original_ppt.columns.to_list(),
             original_ppt.values[0], c='r', alpha=0.45, marker='+',
             label='DWD Obs')  # Obs Same Time (10min difference)

    plt.plot(interpolated_ppt_dwd.columns.to_list(),
             interpolated_ppt_dwd.values[0], c='darkblue', alpha=0.42, marker='*',
             label='DWD Interp RMSE %0.4f' % mse_dwd_interp)

    plt.plot(interpolated_ppt_netatmo_dwd.columns.to_list(),
             interpolated_ppt_netatmo_dwd.values[0], c='darkgreen',
             alpha=0.42, marker='2',
             label='DWD-Netatmo Interp RMSE %0.4f' % mse_dwd_netatmo_interp)

#     plt.plot(interpolated_ppt_netatmo_dwd.columns.to_list(),
#              interpolated_ppt_netatmo_dwd.values[0], c='y',
#              alpha=0.4, marker='3',
#              label='DWD-Netatmo Unc Interp RMSE %0.4f'
#              % mse_dwd_netatmo_unc_interp)

    for col in df_ppt_radolan.columns:
        plt.plot(df_ppt_radolan.index,
                 df_ppt_radolan.loc[:, col].values,
                 c='k', alpha=0.25)

    try:
        plt.plot(df_ppt_radolan.index,
                 df_ppt_radolan.loc[:, col].values,
                 c='k', alpha=0.05,
                 label='Radar RMSE %0.4f' % mse_dwd_radar)
    except Exception as msg:
        print(msg)
        pass
    plt.title('Radolan Station-Adjusted Product: %s \nDWD: %s to %s'
              % (str(radolan_event_date[0]),
                 str(event_start_date[0]),
                 str(event_date_minus_one_hr[0])))
    plt.xticks(rotation=90)
    plt.xlabel('DWD Station ID')
    plt.ylabel('[mm/%s]' % temp_freq)
#     plt.xlabel([df_ppt_obsv.columns.to_list()[::10]], rotation=90)
    plt.legend(loc=0)
    plt.grid(alpha=0.05)
    # plt.tight_layout()
    plt.savefig(
        os.path.join(out_dir,
                     'ppt_event_{}.png'.format(
                         str(event_date[0]).replace(' ',
                                                    '_').replace(':', ''))),
        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)

    plt.close()

    print('Done for this event')
#     break

# rmse per stn
# np.square(np.nanmean(
df_dwd_dwd_rmse_stns = pd.DataFrame.from_dict(
    dwd_dwd_rmse_stns,
    orient='index')

df_dwd_netatmo_rmse_stns = pd.DataFrame.from_dict(
    dwd_netatmo_rmse_stns,
    orient='index')

# df_dwd_netatmo_unc_rmse_stns = pd.DataFrame.from_dict(
#     dwd_netatmo_unc_rmse_stns,
#     orient='index')

df_dwd_radolan_rmse_stns = pd.DataFrame.from_dict(
    dwd_radolan_rmse_stns,
    orient='index')

rmse_stns_dwd_dwd = np.sqrt(
    np.square(df_dwd_dwd_rmse_stns).mean(axis=1)).mean()
rmse_stns_dwd_netatmo = np.sqrt(
    np.square(df_dwd_netatmo_rmse_stns).mean(axis=1)).mean()
# rmse_stns_dwd_netatmo_unc = np.sqrt(
#     np.square(df_dwd_netatmo_unc_rmse_stns).mean(axis=1)).mean()
rmse_stns_dwd_radolan = np.sqrt(
    np.square(df_dwd_radolan_rmse_stns).mean(axis=1)).mean()

df_rmse_all_events = pd.DataFrame(
    index=['each_event_all_stns', 'each_stn_all_events'])

df_rmse_all_events.loc['each_event_all_stns',
                       'dwd_dwd_rmse'] = np.mean(dwd_dwd_rmse)
df_rmse_all_events.loc['each_event_all_stns',
                       'dwd_netatmo_rmse'] = np.mean(dwd_netatmo_rmse)
# df_rmse_all_events.loc['each_event_all_stns',
#                        'dwd_netatmo_unc_rmse'] = np.mean(dwd_netatmo_unc_rmse)
df_rmse_all_events.loc['each_event_all_stns',
                       'dwd_radolan_rmse'] = np.mean(dwd_radolan_rmse)

df_rmse_all_events.loc['each_stn_all_events',
                       'dwd_dwd_rmse'] = rmse_stns_dwd_dwd
df_rmse_all_events.loc['each_stn_all_events',
                       'dwd_netatmo_rmse'] = rmse_stns_dwd_netatmo
# df_rmse_all_events.loc['each_stn_all_events',
#                        'dwd_netatmo_unc_rmse'] = rmse_stns_dwd_netatmo_unc
df_rmse_all_events.loc['each_stn_all_events',
                       'dwd_radolan_rmse'] = rmse_stns_dwd_radolan

df_rmse_all_events.to_csv(
    r'C:\Users\hachem\Desktop\radar\rmse_%s.csv' % temp_freq,
    sep=';')
print('=====\n FINAL RESULTS=====\n')
