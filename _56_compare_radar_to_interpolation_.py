# !/usr/bin/env python.
# -*- coding: utf-8 -*-


# from scipy.interpolate import griddata
from osgeo import osr

import os
# import time
# import timeit

# import wget
# import shapefile
# import numpy as np
import wradlib as wrl
import wradlib.util as util
# import matplotlib.pyplot as plt
# import glob
# import fiona
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
import pandas as pd
# import tables
import numpy as np
import glob
# import datetime
import fiona
import matplotlib.pyplot as plt
import matplotlib
# from matplotlib.patches import Circle, Wedge, Polygon


def find_nearest(array, value):
    ''' given a value, find nearest one to it in original data array'''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


database = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\Final_results\Ppt_ok_ok_un_new2_first_flt__temp_flt__1st_1440min\*.csv'
#base_path = r'C:\Users\hachem\Downloads\SF201907\*'
base_path = r'C:\Users\hachem\Downloads\RW201606\*'

path_to_dwd_ppt = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\all_dwd_hourly_ppt_data_combined_2014_2019_.csv"
path_to_dwd_coords = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\station_coordinates_names_hourly_only_in_BW.csv"

files = glob.glob(base_path)
ppt_data = glob.glob(database)

# hdf5_path = os.path.join(database, 'rain_germany_daily.h5')
# hf = tables.open_file(hdf5_path, 'r')

ishape = fiona.open(
    r"X:\hiwi\ElHachem\Peru_Project\ancahs_dem\DEU_adm\DEU_adm1.shp")
first = ishape.next()

#==============================================================================
#
#==============================================================================
# DWD Coords

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
hourly_events = [  # '2016-06-25 00:00:00',
    '2018-06-11 16:00:00',
    '2018-06-11 17:00:00',
    '2018-06-11 18:00:00',
    '2018-09-23 17:00:00',
    '2018-09-23 18:00:00',
    '2018-09-23 19:00:00']

daily_events = [  # '2018-12-23 00:00:00',
    #'2019-05-22 00:00:00',
    '2018-05-14 00:00:00',
    '2019-07-28 00:00:00']

bound = [0., 1,
         2, 5, 8,
         10, 15, 20,
         25, 30, 35, 40, 45]

cmap = plt.get_cmap('jet')
norm = matplotlib.colors.BoundaryNorm(bound, cmap.N)

for i, file in enumerate(files):
    #     file = file + '.gz'
    print(i, ' / ', len(files), file[50:60])
    event_date = ('20' + file[50:52] + '-' + file[52:54] + '-' + file[54:56]
                  + ' ' + file[56:58] + ':' + file[58:60] + ':00')

    shifted_event = pd.DatetimeIndex([event_date]) + pd.Timedelta(minutes=10)

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

    # convert to lonlat
    radolan_grid_ll = wrl.georef.reproject(radolan_grid_xy,
                                           projection_source=proj_stereo,
                                           projection_target=proj_wgs)
    lon1 = radolan_grid_ll[:, :, 0]
    lat1 = radolan_grid_ll[:, :, 1]

    df_ppt_radolan = pd.DataFrame(index=dwd_in_coords_df.index)
    grid_lons, grid_lats = [], []

    for stn, x0, y0 in zip(dwd_in_coords_df.index,
                           dwd_lons.values, dwd_lats.values):
        #         print(x0, y0)
        grd_lon_loc = find_nearest(lon1.flatten(), x0)
        grd_lat_loc = find_nearest(lat1.flatten(), y0)
        grid_lons.append(grd_lon_loc)
        grid_lats.append(grd_lat_loc)
#         plt.ioff()
#         plt.scatter(x0, y0)
#         plt.scatter(grd_lon_loc, grd_lat_loc)
#         plt.show()
        ix_lon_loc = np.where(lon1.flatten() == grd_lon_loc)[0]
        ix_lat_loc = np.where(lat1.flatten() == grd_lat_loc)[0]

        ix_lon_in_arr = np.unravel_index(ix_lon_loc, (900, 900))
#         lon1[]
        ix_lat_in_arr = np.unravel_index(ix_lat_loc, (900, 900))

        radolan_ppt_loc0 = rwdata.data[ix_lon_in_arr[0], ix_lat_in_arr[0]]
        radolan_ppt_loc1 = rwdata.data[ix_lon_in_arr[1], ix_lat_in_arr[1]]

        df_ppt_radolan.loc[stn, 'loc0'] = radolan_ppt_loc0
        df_ppt_radolan.loc[stn, 'loc1'] = radolan_ppt_loc1


#         print(ix_lon_in_arr[1], ix_lat_in_arr[1])

    df_ppt_radolan[df_ppt_radolan < 0] = np.nan
    df_ppt_radolan.dropna(how='all')

    # obsv
    df_ppt_obsv = dwd_in_ppt_vals_df.loc[shifted_event, :]

    plt.ioff()

    plt.figure()
    plt.plot(df_ppt_obsv.columns.to_list(),
             df_ppt_obsv.values[0], c='r', alpha=0.75,
             label='Obs')
    plt.plot(df_ppt_radolan.index, df_ppt_radolan.loc0.values, c='b', alpha=0.75,
             label='Radar loc0')
    plt.plot(df_ppt_radolan.index, df_ppt_radolan.loc1.values, c='g', alpha=0.75,
             label='Radar loc1')
    plt.xticks([])
    plt.ylabel('mm/h')
#     plt.xlabel([df_ppt_obsv.columns.to_list()[::10]], rotation=90)
    plt.legend(loc=0)
    plt.tight_layout()
    plt.savefig(
        os.path.join(r'X:\exchange\ElHachem\Figures_Netatmo\new_events\radar',
                     'ppt_event_{}.png'.format(file[49:59])), dpi=600)

#     plt.show()
#     break
    plt.close()


#         radolan_ppt.append(radolan_ppt_loc1)

    mask = np.ones_like(lon1, dtype=np.bool)

    for n, i_poly in enumerate(first['geometry']['coordinates']):
        # print(n)
        #         plt.plot(np.array(i_poly)[0, :, 0],
        #                  np.array(i_poly)[0, :, 1],
        #                  linestyle='-',
        #                  linewidth=0.8,
        #                  color='black'
        #                  )
        p = path.Path(np.array(i_poly)[0, :, :])
        grid_mask = p.contains_points(
            np.vstack((lon1.flatten(), lat1.flatten())).T).reshape(900, 900)
        mask[grid_mask] = 0

    rwdata[mask] = -1
    rw_maskes = np.ma.masked_array(rwdata, rwdata < 0.)

    plt.figure()
#     min_x = xss[np.argmin([xs - x0 for xs in xss])]
#     min_y = yss[np.argmin([ys - y0 for ys in yss])]

    plt.pcolormesh(lon1, lat1, rw_maskes, cmap=cmap,
                   vmin=0, norm=norm)

    plt.ylabel('Latitude [°]')
    plt.xlabel('Longitude [°]')

    plt.xlim([7.1, 10.7])
    plt.ylim([47.3, 50.0])


#     plt.axis('equal')
    # radar.set_aspect('equal', 'box')
    # plt.xlim([netatmo.get_xlim()[0], netatmo.get_xlim()[1]])
    # plt.ylim([netatmo.get_ylim()[0], netatmo.get_ylim()[1]])
    cbar = plt.colorbar()
    # cbar.set_label('Hourly Precipitation [mm]', rotation=270)
    # colorbar.colorbar(ax0)
#     plt.axis('equal')
    plt.scatter(dwd_lons.values, dwd_lats.values,
                marker='x', alpha=0.5, color='k', s=10)
    plt.scatter(grid_lons, grid_lats, marker='.', color='r', alpha=0.5, s=10)

    plt.tight_layout()
    plt.savefig(
        os.path.join(r'X:\exchange\ElHachem\Figures_Netatmo\new_events\radar',
                     'event_{}.png'.format(file[49:59])), dpi=600)

#     plt.show()
#     break
    plt.close()
    # plt.show()
