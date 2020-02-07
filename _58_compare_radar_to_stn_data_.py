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

# # from matplotlib.patches import Circle, Wedge, Polygon
#
#


def find_nearest_nbrs(array, value, nbr_ngbrs=9):
    ''' given a value, find nearest one to it in original data array'''
    array = np.asarray(array)
    idx = np.argsort((np.abs(array - value)))[:nbr_ngbrs]  # .argmin()
    return array[idx]


# mas for BW data
mask = np.load(r'C:\Users\hachem\Desktop\radar\masked_array_bw',
               allow_pickle=True)

#
#
database = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
            r'\oridinary_kriging_compare_DWD_Netatmo\Final_results'
            r'\Ppt_ok_ok_un_new2_first_flt__temp_flt__1st_1440min\*.csv')

base_path = r'C:\Users\hachem\Downloads\SF201907\*'
#base_path = r'C:\Users\hachem\Downloads\SF201907\*'
# base_path = r'C:\Users\hachem\Downloads\RW201809\*'

path_to_dwd_ppt = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                   r"\NetAtmo_BW\ppt_all_dwd_1440min_.csv")
path_to_dwd_coords = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                      r"\NetAtmo_BW\station_coordinates_names_hourly_only_in_BW.csv")
#
files = glob.glob(base_path)
ppt_data = glob.glob(database)
#


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

event_date_initial = ('20' + files[0][50:52] + '-' + files[0][52:54] +
                      '-' + files[0][54:56] + ' ' + files[0][56:58] + ':' +
                      files[0][58:60] + ':00')

for i, file in enumerate(files):
    #     file = file + '.gz'

    event_date = ('20' + file[50:52] + '-' + file[52:54] + '-' + file[54:56]
                  + ' ' + file[56:58] + ':' + file[58:60] + ':00')

    print('Event: ', i, '/', len(files),  event_date)

    event_date = pd.DatetimeIndex([event_date])

    if event_date.hour == 23:
        #         if event_date in hourly_events:
        event_date_minus_one_hr = event_date - pd.Timedelta(minutes=60)

        event_time = event_date + pd.Timedelta(minutes=10)
        shifted_event = event_date_minus_one_hr + pd.Timedelta(minutes=10)

#         shifted_event_str = str(shifted_event[0]).replace('23', '00')
#         dwd_in_ppt_vals_df.loc[shifted_event_str,:]
    #     shifted_event = event_time
        rwdata, rwattrs = wrl.io.read_radolan_composite(file)
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

        df_ppt_radolan = pd.DataFrame(index=dwd_in_coords_df.index)

        ix_lon_loc = []
        ix_lat_loc = []
        texts = []

        plt.ioff()
        fig = plt.figure(figsize=(20, 12), dpi=200)
        ax = fig.add_subplot(111)

        print('\nGetting station data\n')

        for stn, x0, y0 in zip(dwd_in_coords_df.index,
                               dwd_lons.values, dwd_lats.values):

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
        df_ppt_obsv_orig = dwd_in_ppt_vals_df.loc[
            event_time, :].dropna(how='all')
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
            Path(os.path.join(
                r'C:\Users\hachem\Desktop\radar\ppt_radar_vs_ppt_dwd',
                'radar_event_{}.png'.format(file[50:59]))))

    #     plt.show()

        plt.close()
        print('Plotting station data\n')
        #======================================================================
        #
        #======================================================================
        plt.figure(figsize=(20, 12), dpi=200)

        plt.plot(df_ppt_obsv_orig.columns.to_list(),
                 df_ppt_obsv_orig.values[0], c='g', alpha=0.75,
                 label='DWD stns')  # Obs Same Time (10min difference)

#         plt.plot(df_ppt_obsv_shifted.columns.to_list(),
#                  df_ppt_obsv_shifted.values[0], c='b', alpha=0.75,
#                  label='Obs Shifted Time (-1hr)')

        for col in df_ppt_radolan.columns:
            plt.plot(df_ppt_radolan.index,
                     df_ppt_radolan.loc[:, col].values,
                     c='k', alpha=0.25)

        plt.plot(df_ppt_radolan.index,
                 df_ppt_radolan.loc[:, col].values,
                 c='k', alpha=0.25, label='Radar')
        plt.title(str(event_date[0]))
        plt.xticks(rotation=90)
        plt.ylabel('[mm/day]')
    #     plt.xlabel([df_ppt_obsv.columns.to_list()[::10]], rotation=90)
        plt.legend(loc=0)
        # plt.tight_layout()
        plt.savefig(
            os.path.join(r'C:\Users\hachem\Desktop\radar\ppt_radar_vs_ppt_dwd',
                         'ppt_event_{}.png'.format(file[50:59])))

        plt.close()

        print('Done for this event')
    #     break
