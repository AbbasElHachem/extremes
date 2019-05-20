'''
Created on %(date)s

@author: %(username)s
'''

import os
import timeit
import time
import pandas as pd
import fiona
import numpy as np
import shapely.geometry as shg
import utm


from netatmo_core import (
    get_auth,
    get_public_data,
    get_measured_data,
    assemble_stn_measure_data,
    refresh_token)


def main():

    main_dir = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW'
    os.chdir(main_dir)

    # os.chdir(main_dir)

    payload_auth = {
        'grant_type': 'password',
        # 'username': 'micha.eisele@iws.uni-stuttgart.de',
        # 'password': 'QWEasd123!',
        #         'username': 'faizananwar2006@gmail.com',
        #         'password': 'n37AT30+B',
        'username': 'abbas.el-hachem@iws.uni-stuttgart.de',
        'password': 'Netatmo159',
        'client_id': '59788efd6b0affa5208b480e',
        'client_secret': 'ZtvUyd2q32YIJrV4EC1jGyatm',
        'scope': 'read_station'}

    gauge_type = 'rain'
#     gauge_type = ' humidity'
    # BW
    public_data_params = {
        'lat_ne': '49.847',
        'lon_ne': '10.501',
        'lat_sw': '47.571',
        'lon_sw': '7.494',
        'required_data': gauge_type}

    # load shapefile
    location = 'bw'
    path_shapefile = r'X:\hiwi\ElHachem\GitHub\extremes\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89.shp'

#     # Germany
#     public_data_params = {
#         'lat_ne': '55.503750',
#         'lon_ne': '14.633789',
#         'lat_sw': '47.070122',
#         'lon_sw': '5.361328',
#         'required_data': gauge_type}
#
#     # load shapefile
#     location = 'germany'
#     path_shapefile = r'P:\2017_SYNOPSE_II\02_Import\03_QGIS_project_germany\germany_shape_utm32.shp'

    beg_date = '2009-01-01'
    end_date = '2019-03-31'
    in_date_fmt = '%Y-%m-%d'

    # temporal scale
    temp_scale = '1hour'
    # temp_scale = '5min'
#     temp_scale = '1day'

    scale_type = 'sum_rain'
#     scale_type = 'max_hum'

    data_type = 'rain_60min'
#     data_type = 'date_max_hum'

    overwrite_flag = False

    ishape = fiona.open(path_shapefile)
    first = ishape.next()
    polygon = np.array(first['geometry']['coordinates'][0])
    try:
        bawue = polygon[0, :, :]
    except IndexError:
        bawue = polygon[:, :]

    if location == 'germany':
        bawue = np.array(first['geometry']['coordinates'][-1])[0, :, :]

    poly = shg.Polygon(bawue)

    out_dir_daily = '%s_%s_%s' % (gauge_type, location, temp_scale)

    out_daily_data_loc = os.path.join(
        out_dir_daily, 'netatmo_%s_%s_((%s)).csv')

    out_coords_loc = os.path.join(
        out_dir_daily, 'netatmo_%s_%s_coords.csv' % (location, temp_scale))

    if not os.path.exists(out_dir_daily):
        os.makedirs(out_dir_daily)

    # write header to coordinates file
    if os.path.isfile(out_coords_loc) is False:
        f = open(out_coords_loc, 'w')
        f.write(str('stn_mac; lon; lat; elev; module_mac; tz'))
        f.write('\n')
        f.close()

    sep = ';'

    beg_date_ts = int(time.mktime(time.strptime(beg_date, in_date_fmt)))
    end_date_ts = int(time.mktime(time.strptime(end_date, in_date_fmt)))

    print('Authorizing...')
    response_auth = get_auth(payload_auth)
    print(response_auth)

    public_data_params['access_token'] = response_auth['access_token']

    print('Getting public data...')
    public_data = get_public_data(
        public_data_params, payload_auth, response_auth)

    # device_module_id_dict = {}
    data_daily_time_ser_dict = {}
    # device_id_meta_data_dict = {}

    print('Found %d stns!' % len(public_data))

    if not public_data:
        print('No public data!')
        return

    print('\n' * 3)
    print('Getting measured data...')

    for stn_no, stn_meta_data in enumerate(public_data):
        # check if station within shapefile
        utm_val = utm.from_latlon(stn_meta_data['place']['location'][1],
                                  stn_meta_data['place']['location'][0])
        if poly.contains(shg.Point(tuple(utm_val[:2]))) == True:

            stn_mac_add = stn_meta_data['_id']

            _out_loc = out_daily_data_loc % (location, temp_scale, stn_mac_add)
            _out_loc = _out_loc.replace(':', '_')

            print('Stn_no:', stn_no + 1)
            print(stn_meta_data)

            if os.path.exists(_out_loc) and (not overwrite_flag):
                print('File: %s exists. Not overwriting!' % _out_loc)
                continue

            # device_id_meta_data_dict[stn_mac_add] = stn_meta_data

            measures_dict = stn_meta_data['measures']
            module_mac = ''
            for _ in measures_dict:
                _keys = list(measures_dict[_].keys())
                if 'type' in _keys:
                    if data_type in list(measures_dict[_]['type']):
                        module_mac = _
                        break

                elif data_type in list(measures_dict[_].keys()):
                    module_mac = _
                    break

            if not module_mac:
                print(
                    ('Stn %s did not have any valid module with the given '
                     'data_type!') % stn_mac_add)

                print('The returned measure_dict is:', measures_dict, '\n')
                continue

            curr_beg_date_ts = beg_date_ts
            curr_data_daily_time_ser_list = []

            while curr_beg_date_ts < end_date_ts:

                measure_data_params = {
                    'access_token': response_auth['access_token'],
                    'device_id': stn_mac_add,
                    'module_id': module_mac,
                    'scale': temp_scale,
                    'type': scale_type,
                    'real_time': 'true',
                    'date_begin': str(curr_beg_date_ts),
                    'dashboard_data': 'Rain'}

                refresh_token(payload_auth, response_auth)

                measure_data = get_measured_data(
                    measure_data_params, payload_auth, response_auth)

                # device_module_id_dict[stn_mac_add] = module_mac

                if measure_data:
                    curr_data_daily_time_ser_list.append(
                        assemble_stn_measure_data(measure_data, stn_mac_add))

                    curr_end_date = curr_data_daily_time_ser_list[-1].index[-1]

                    print(
                        'Total steps: %d' %
                        curr_data_daily_time_ser_list[-1].shape[0])

                    curr_beg_date_ts = int(curr_end_date.timestamp())

                else:
                    print('\nCould not get more data for this station!')

                    print(
                        'End date was:',
                        pd.to_datetime(curr_beg_date_ts, unit='s'))

                    print('measure_data:', measure_data)
                    break

            if curr_data_daily_time_ser_list:
                curr_stn_daily_comb_ser = pd.concat(
                    curr_data_daily_time_ser_list)

                curr_stn_daily_comb_ser.sort_index(inplace=True)

                curr_stn_daily_comb_ser.to_csv(
                    # _out_loc, sep=sep, float_format='%0.2f')
                    _out_loc, sep=sep, float_format='%0.5f')

                data_daily_time_ser_dict[stn_mac_add] = curr_stn_daily_comb_ser

                # add meta data to coord file
                if 'place' in stn_meta_data:
                    curr_coords_data = stn_meta_data['place']

                    if 'location' in curr_coords_data:
                        curr_lon, curr_lat = curr_coords_data['location']
                    else:
                        curr_lon, curr_lat = [pd.np.nan] * 2

                    if 'altitude' in curr_coords_data:
                        curr_elev = curr_coords_data['altitude']
                    else:
                        curr_elev = pd.np.nan

                    if 'timezone' in curr_coords_data:
                        curr_tz = curr_coords_data['timezone']
                    else:
                        curr_tz = pd.np.nan

                else:
                    [curr_lon,
                     curr_lat,
                     curr_elev,
                     curr_tz] = [pd.np.nan] * 4

                # [stn_mac_add, curr_lon, curr_lat, curr_elev, module_mac, curr_tz]
                # ['float', 'float', 'int', 'object', 'object']

                f = open(out_coords_loc, 'a')
                f.write(str('%s; %0.8f; %0.8f; %i; %s; %s') % (stn_mac_add, curr_lon, curr_lat,
                                                               curr_elev, module_mac, curr_tz))
                f.write('\n')
                f.close()

                print('Wait 5 secs...')
                time.sleep(5)

            print('\n')

    print('\n' * 3)

    return


if __name__ == '__main__':
    print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
    # START = timeit.default_timer()  # to get the runtime of the program

    main()

    # STOP = timeit.default_timer()  # Ending time
    # print(('\n\a\a\a Done with everything on %s. Total run time was'
    #        ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
