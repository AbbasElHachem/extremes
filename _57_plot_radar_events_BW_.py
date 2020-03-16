# !/usr/bin/env python.
# -*- coding: utf-8 -*-


import os

import glob
import fiona
import matplotlib
import pandas as pd
import numpy as np
import wradlib as wrl
import matplotlib.pyplot as plt

from matplotlib import path
from osgeo import osr

from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap


#base_path = r'C:\Users\hachem\Downloads\SF201907\*'
base_path = r'C:\Users\hachem\Downloads\SF201407.tar\SF201407\*'
state_name = 'Rheinland-Pfalz'

files = glob.glob(base_path)

# hdf5_path = os.path.join(database, 'rain_germany_daily.h5')
# hf = tables.open_file(hdf5_path, 'r')

ishape = (r"P:\2017_SYNOPSE_II\02_Import\07_Shapefiles\DEU_adm_shp"
          r"\DEU_adm1.shp")


#first = ishape.next()


shp_objects_all = [shp for shp in list(fiona.open(ishape))
                   if shp['properties']['NAME_1'] == state_name]

#==============================================================================
#
# #============================================================================
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

bound = [0., 1,
         2, 5, 8,
         10, 15, 20,
         25, 30]  # , 35, 40, 45]
interval_ppt = np.linspace(0.05, 0.95)
colors_ppt = plt.get_cmap('jet_r')(interval_ppt)
cmap_ppt = LinearSegmentedColormap.from_list('name', colors_ppt)
#cmap_ppt = plt.get_cmap('jet_r')
cmap_ppt.set_over('navy')

norm = matplotlib.colors.BoundaryNorm(bound, cmap_ppt.N)

for i, file in enumerate(files):
    #     file = file + '.gz'

    #     event_date = ('20' + file[50:52] + '-' + file[52:54] + '-' + file[54:56]
    #                   + ' ' + file[56:58] + ':' + file[58:60] + ':00')
    #     print(i, '/', len(files),  event_date)
    # #     if event_date in hourly_events:
    #     event_date_minus_one_hr = pd.DatetimeIndex(
    #         [event_date]) - pd.Timedelta(minutes=60)

    #     shifted_event = event_date_minus_one_hr + pd.Timedelta(minutes=10)

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

    mask = np.ones_like(lon1, dtype=np.bool)
    # first['geometry']['coordinates']
    for n, i_poly_all in enumerate(shp_objects_all):
        i_poly = i_poly_all['geometry']['coordinates']
        if 0 < len(i_poly) <= 1:
            p = path.Path(np.array(i_poly)[0])
            grid_mask = p.contains_points(
                np.vstack((lon1.flatten(),
                           lat1.flatten())).T).reshape(900, 900)
            mask[grid_mask] = 0
        else:
            for ix in range(len(i_poly)):

                p = path.Path(np.array(i_poly[ix]))
                grid_mask = p.contains_points(
                    np.vstack((lon1.flatten(),
                               lat1.flatten())).T).reshape(900, 900)
                mask[grid_mask] = 0
#     mask.dump(r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz\mask.npy')
#     rwdata[mask] = -1
    rw_maskes = np.ma.masked_array(rwdata, rwdata < 0.)

#         plt.figure()
#     #     min_x = xss[np.argmin([xs - x0 for xs in xss])]
#     #     min_y = yss[np.argmin([ys - y0 for ys in yss])]
#
#         plt.pcolormesh(lon1, lat1, rw_maskes, cmap=cmap_ppt,
#                        vmin=0, norm=norm, vmax=30)
#
#         plt.ylabel('Latitude [°]')
#         plt.xlabel('Longitude [°]')
#
#         plt.xlim([7.1, 10.7])
#         plt.ylim([47.3, 50.0])
#
#     #     plt.axis('equal')
#         # radar.set_aspect('equal', 'box')
#         # plt.xlim([netatmo.get_xlim()[0], netatmo.get_xlim()[1]])
#         # plt.ylim([netatmo.get_ylim()[0], netatmo.get_ylim()[1]])
#         cbar = plt.colorbar(extend='max', label='[mm/h]')
#
#         plt.tight_layout()
#         plt.savefig(
#             Path(os.path.join(
#                 r'C:\Users\hachem\Desktop\radar',
#                 'hourly_event_{}.png'.format(file[50:59]))), dpi=600)
#
#     #     plt.show()
#     #     break
#         plt.close()
