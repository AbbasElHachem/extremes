# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
Construct DISTANCE Matrix between Netamo Rainfall and 
    Netamo Humidity Stations based on coordinates 
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


# coords_netatmo_humidity_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#                                    r'\NetAtmo_BW\humidity_bw_1hour'
#                                    r'\netatmo_bw_1hour_coords.csv')
# assert os.path.exists(coords_netatmo_humidity_df_file), 'wrong Hum coords file'

# coords_netatmo_temp_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
#                                r'\NetAtmo_BW\temperature_bw_1hour'
#                                r'\netatmo_bw_1hour_coords.csv')
# assert os.path.exists(coords_netatmo_temp_df_file), 'wrong Temp coords file'

coords_netatmo_temp_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                               r'\NetAtmo_BW\rain_bw_1hour'
                               r'\netatmo_bw_1hour_coords.csv')

assert os.path.exists(coords_netatmo_temp_df_file), 'wrong NETATMO coords file'

coords_netatmo_ppt_df_file = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                              r'\NetAtmo_BW\rain_bw_1hour'
                              r'\netatmo_bw_1hour_coords.csv')

assert os.path.exists(coords_netatmo_ppt_df_file), 'wrong NETATMO coords file'

out_save_dir = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW')


if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)

lon_col_name = 'lon'
lat_col_name = 'lat'

df_sep = ';'

# def epsg wgs84 and utm32 for coordinates conversion
wgs82 = "+init=EPSG:4326"
utm32 = "+init=EPSG:32632"

in_df_ppt_netatmo_coords = pd.read_csv(coords_netatmo_ppt_df_file,
                                       index_col=0,
                                       sep=df_sep)

x_ppt_netatmo, y_ppt_netatmo = convert_coords_fr_wgs84_to_utm32_(
    wgs82, utm32,
    in_df_ppt_netatmo_coords.loc[:, lon_col_name].values.ravel(),
    in_df_ppt_netatmo_coords.loc[:, lat_col_name].values.ravel())
#==============================================================================

in_df_temp_netatmo_coords = pd.read_csv(coords_netatmo_temp_df_file,
                                        index_col=0,
                                        sep=df_sep)
in_df_temp_netatmo_coords.drop_duplicates(keep='first', inplace=True)

x_hum_netatmo, y_hum_netatmo = convert_coords_fr_wgs84_to_utm32_(
    wgs82, utm32,
    in_df_temp_netatmo_coords.loc[:, lon_col_name].values.ravel(),
    in_df_temp_netatmo_coords.loc[:, lat_col_name].values.ravel())
#==============================================================================


coords_ppt_netatmo = np.array([(x, y) for x, y in zip(x_ppt_netatmo,
                                                      y_ppt_netatmo)])
coords_hum_netatmo = np.array([(x, y) for x, y in zip(x_hum_netatmo,
                                                      y_hum_netatmo)])


# distance.cdist(coords_dwd, coords_netatmo, 'euclidean')

data_mtx = np.zeros(
    shape=(in_df_ppt_netatmo_coords.index.shape[0],
           in_df_temp_netatmo_coords.index.shape[0])).astype('float')
data_mtx[data_mtx == 0] = np.nan

stations_ppt_id_with_ = [stn_id.replace(':', '_')
                         for stn_id in in_df_ppt_netatmo_coords.index]

stations_hum_id_with_ = [stn_id.replace(':', '_')
                         for stn_id in in_df_temp_netatmo_coords.index]
df_distance = pd.DataFrame(index=stations_ppt_id_with_,
                           columns=stations_hum_id_with_,
                           data=data_mtx)

for stn_mac in in_df_ppt_netatmo_coords.index:
    try:
        print(stn_mac)
        lon_stn = in_df_ppt_netatmo_coords.loc[stn_mac, lon_col_name]
        lat_stn = in_df_ppt_netatmo_coords.loc[stn_mac, lat_col_name]
        assert isinstance(lon_stn, np.float)
    except Exception as msg:
        print(stn_mac, msg, '\n Coordinates not unique')
        lon_stn = lon_stn.values[0]
        lat_stn = lat_stn.values[0]

    try:
        xstn, ystn = convert_coords_fr_wgs84_to_utm32_(wgs82,
                                                       utm32,
                                                       lon_stn,
                                                       lat_stn)

        dist = distance.cdist([(xstn, ystn)], coords_hum_netatmo, 'euclidean')
        stn_max_with_ = stn_mac.replace(':', '_')
        df_distance.loc[stn_max_with_, :] = dist
    except Exception as msg:
        print(msg)
        raise Exception

df_distance.to_csv(os.path.join(out_save_dir,
                                'distance_mtx_in_m_NetAtmo_ppt_Netatmo_ppt.csv'),
                   sep=';')

print('Done with everything')
