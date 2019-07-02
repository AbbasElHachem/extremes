# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Construct DISTANCE Matrix between Netamo and DWD Stations based on coordinates 
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"
#==============================================================================
#
#==============================================================================
import os

import numpy as np
import pandas as pd


from scipy.spatial import distance
from _00_additional_functions import (convert_coords_fr_wgs84_to_utm32_)


#==============================================================================
#
#==============================================================================


coords_dwd_df_file = (r'X:\exchange\ElHachem\niederschlag_deutschland'
                      r'\1993_2016_5min_merge_nan.csv')
assert os.path.exists(coords_dwd_df_file), 'wrong DWD coords file'

# coords_netatmo_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#                           r'\NetAtmo_BW\rain_bw_1hour\netatmo_bw_1hour_coords.csv')
#
# assert os.path.exists(coords_netatmo_df_file), 'wrong NETATMO coords file'

# coords_netatmo_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#                           r'\NetAtmo_BW\rain_bw_5min\netatmo_bw_1hour_coords.csv')

coords_netatmo_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                          r'\NetAtmo_BW\\temperature_bw_1hour\netatmo_bw_1hour_coords.csv')

assert os.path.exists(coords_netatmo_df_file), 'wrong NETATMO coords file'

out_save_dir = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW')


if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)

lon_col_name = ' lon'
lat_col_name = ' lat'

x_col_name = 'Rechtswert'
y_col_name = 'Hochwert'

# def epsg wgs84 and utm32 for coordinates conversion
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

in_df_netatmo_coords = pd.read_csv(coords_netatmo_df_file, index_col=0,
                                   sep=';')
# drop all duplicates in stations
in_df_netatmo_coords.drop_duplicates(keep='first', inplace=True)

xnetatmo, ynetatmo = convert_coords_fr_wgs84_to_utm32_(
    wgs82, utm32,
    in_df_netatmo_coords.loc[:, lon_col_name].values.ravel(),
    in_df_netatmo_coords.loc[:, lat_col_name].values.ravel())

in_df_dwd_coords = pd.read_csv(coords_dwd_df_file, index_col=3, sep=';')

x_dwd, y_dwd = (in_df_dwd_coords.loc[:, x_col_name].values.ravel(),
                in_df_dwd_coords.loc[:, y_col_name].values.ravel())


coords_netatmo = np.array([(x, y) for x, y in zip(xnetatmo, ynetatmo)])
coords_dwd = np.array([(xd, yd) for xd, yd in zip(x_dwd, y_dwd)])


# distance.cdist(coords_dwd, coords_netatmo, 'euclidean')

data_mtx = np.zeros(
    shape=(in_df_netatmo_coords.index.shape[0],
           in_df_dwd_coords.index.shape[0])).astype('float')
data_mtx[data_mtx == 0] = np.nan

stations_id_with_ = [stn_id.replace(':', '_')
                     for stn_id in in_df_netatmo_coords.index]

df_distance = pd.DataFrame(index=stations_id_with_,
                           columns=in_df_dwd_coords.index,
                           data=data_mtx)

for stn_mac in in_df_netatmo_coords.index:
    try:
        print(stn_mac)
        lon_stn = in_df_netatmo_coords.loc[stn_mac, lon_col_name]
        lat_stn = in_df_netatmo_coords.loc[stn_mac, lat_col_name]
        xstn, ystn = convert_coords_fr_wgs84_to_utm32_(wgs82,
                                                       utm32,
                                                       lon_stn,
                                                       lat_stn)

        dist = distance.cdist([(xstn, ystn)], coords_dwd, 'euclidean')
        stn_max_with_ = stn_mac.replace(':', '_')
        df_distance.loc[stn_max_with_, :] = dist
        print('done for stn', stn_mac)
    except Exception as msg:
        print(msg, 'error for', stn_mac)
        continue

df_distance.to_csv(os.path.join(out_save_dir,
                                'distance_mtx_in_m_NetAtmo_temp_DWD_ppt.csv'),
                   sep=';')

print('Done with everything')
