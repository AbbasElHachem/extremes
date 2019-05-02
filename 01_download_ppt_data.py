'''
Created on %(date)s

@author: %(username)s
'''

import os
import timeit
import time

import pandas as pd

from netatmo_core import (
    get_auth,
    get_public_data,
    get_measured_data,
    assemble_stn_measure_data,
    refresh_token)


def main(idx):

    #     main_dir = os.path.join(os.path.dirname(__file__))
    main_dir = r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW'
    os.chdir(main_dir)

    payload_auth = {
        'grant_type': 'password',
        'username': 'faizananwar2006@gmail.com',
        'password': 'n37AT30+B',
        'client_id': '59788efd6b0affa5208b480e',
        'client_secret': 'ZtvUyd2q32YIJrV4EC1jGyatm',
        'scope': 'read_station'}

    gage_type = 'rain'

    public_data_params = {
        'lat_ne': '49.847',
        'lon_ne': '10.501',
        'lat_sw': '47.571',
        'lon_sw': '7.494',
        'required_data': gage_type}

    beg_date = '2009-01-01'
    end_date = '2019-03-31'
    in_date_fmt = '%Y-%m-%d'

    scale = '1hour'
    data_type = 'rain_60min'

    down_events_flag = False
    overwrite_flag = False

    min_ppt_thresh = 1.0

    out_dir_daily = 'ppt_bw_grosser_daily'
    out_dir_sub_daily = 'ppt_bw_grosser_sub_daily'

    out_daily_data_loc = os.path.join(
        out_dir_daily, 'netatmo_bw_ppt__((%s)).csv')

    out_sub_daily_data_loc = os.path.join(
        out_dir_sub_daily, 'netatmo_bw_ppt__((%s)).csv')

    out_coords_loc = os.path.join(
        out_dir_daily, 'netatmo_bw_ppt_coords_%d.csv' % idx)

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

    device_module_id_dict = {}
    data_daily_time_ser_dict = {}
    device_id_meta_data_dict = {}

    if down_events_flag:
        data_sub_daily_time_ser_dict = {}

    print('Found %d stns!' % len(public_data))

    if not public_data:
        print('No public data!')
        return

    if not os.path.exists(out_dir_daily):
        os.makedirs(out_dir_daily)

    if (not os.path.exists(out_dir_sub_daily)) and down_events_flag:
        os.makedirs(out_dir_sub_daily)

    print('\n' * 3)
    print('Getting measured data...')

    for stn_no, stn_meta_data in enumerate(public_data):
        stn_mac_add = stn_meta_data['_id']

        _out_loc = out_daily_data_loc % stn_mac_add
        _out_loc = _out_loc.replace(':', '_')

        print('Stn_no:', stn_no + 1)
        print(stn_meta_data)

        if os.path.exists(_out_loc) and (not overwrite_flag):
            print('File: %s exists. Not overwriting!' % _out_loc)
            continue

        device_id_meta_data_dict[stn_mac_add] = stn_meta_data

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
                'scale': scale,
                'type': gage_type,
                'real_time': 'true',
                'date_begin': str(curr_beg_date_ts)}

            refresh_token(payload_auth, response_auth)

            measure_data = get_measured_data(
                measure_data_params, payload_auth, response_auth)

            device_module_id_dict[stn_mac_add] = module_mac

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
            curr_stn_daily_comb_ser = pd.concat(curr_data_daily_time_ser_list)

            curr_stn_daily_comb_ser.sort_index(inplace=True)

            curr_stn_daily_comb_ser.to_csv(
                _out_loc, sep=sep, float_format='%0.2f')

            data_daily_time_ser_dict[stn_mac_add] = curr_stn_daily_comb_ser
            print('Wait 5 secs...')
            time.sleep(5)

        print('\n' * 3)

    all_stn_macs = list(device_module_id_dict.keys())

    coords_cols = ['lon', 'lat', 'elev', 'module_mac', 'tz']
    coords_dtypes = ['float', 'float', 'int', 'object', 'object']
    cols_dtypes_dict = dict(zip(coords_cols, coords_dtypes))

    all_stns_coords_df = pd.DataFrame(
        index=all_stn_macs, columns=coords_cols, dtype=object)

    for _stn in all_stn_macs:
        curr_meta_data = device_id_meta_data_dict[_stn]
        curr_mod_mac = device_module_id_dict[_stn]

        if 'place' in curr_meta_data:
            curr_coords_data = curr_meta_data['place']

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

        all_stns_coords_df.loc[_stn, :] = [
            curr_lon, curr_lat, curr_elev, curr_mod_mac, curr_tz]

    all_stns_coords_df = all_stns_coords_df.astype(dtype=cols_dtypes_dict)
    all_stns_coords_df.to_csv(
        out_coords_loc, sep=sep, float_format='%0.8f', index_label='stn_mac')

    print('\n' * 3)
    if down_events_flag:
        print('Downloading sub daily ppt event data...')

        for curr_stn_mac_add in data_daily_time_ser_dict:
            print('Current station:', curr_stn_mac_add)

            curr_module_mac = device_module_id_dict[curr_stn_mac_add]

            curr_stn_sub_daily_time_ser_list = []

            curr_daily_ser = data_daily_time_ser_dict[curr_stn_mac_add]

            for j, date in enumerate(curr_daily_ser.index):
                if not j:
                    continue

                if curr_daily_ser.loc[
                        date, curr_stn_mac_add] <= min_ppt_thresh:

                    continue

                curr_beg_date_ts = int(
                    curr_daily_ser.index[j - 1].timestamp())

                curr_end_date_ts = int(curr_daily_ser.index[j].timestamp())

                measure_data_params = {
                    'access_token': response_auth['access_token'],
                    'device_id': curr_stn_mac_add,
                    'module_id': curr_module_mac,
                    'scale': 'max',
                    'type': 'rain',
                    'real_time': 'true',
                    'date_begin': str(curr_beg_date_ts),
                    'date_end': str(curr_end_date_ts)}

                measure_data = get_measured_data(
                    measure_data_params, payload_auth, response_auth)

                if measure_data:
                    curr_stn_sub_daily_time_ser_list.append(
                        assemble_stn_measure_data(
                            measure_data, curr_stn_mac_add))

                    print(
                        'Total steps: %d' %
                        curr_stn_sub_daily_time_ser_list[-1].shape[0])

                else:
                    print('\nCould not get more data for this station!')
                    print('measure_data:', measure_data)
                    break

            if curr_stn_sub_daily_time_ser_list:
                curr_stn_sub_daily_comb_ser = pd.concat(
                    curr_stn_sub_daily_time_ser_list)

                _out_loc = out_sub_daily_data_loc % curr_stn_mac_add
                _out_loc = _out_loc.replace(':', '_')

                curr_stn_sub_daily_comb_ser.sort_index(inplace=True)

                curr_stn_sub_daily_comb_ser.to_csv(
                    _out_loc, sep=sep, float_format='%0.2f')

                data_sub_daily_time_ser_dict[curr_stn_mac_add] = (
                    curr_stn_sub_daily_comb_ser)

                print('Wait 5 secs...')
                time.sleep(5)

            print('\n' * 3)

    return


if __name__ == '__main__':
    print('\a\a\a\a Started on %s \a\a\a\a\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program
#     for idx in range(2):
    main(1)

    STOP = timeit.default_timer()  # Ending time
    print(('\n\a\a\a Done with everything on %s. Total run time was'
           ' about %0.4f seconds \a\a\a' % (time.asctime(), STOP - START)))
