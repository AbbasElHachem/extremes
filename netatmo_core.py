"""
Created on Wed Aug  2 09:02:20 2017

@author: Faizan
"""

import timeit
import time
import requests

import pandas as pd

LAST_REFRESH_AUTH_AT = None
AUTH_EXP_IN = None

TEN_SEC_REQ_CTR = 0
LAST_TEN_SEC_REQ_AT = timeit.default_timer()
CURR_TEN_SEC_REQ_AT = timeit.default_timer()

ONE_HOUR_REQ_CTR = 0
LAST_ONE_HOUR_REQ_AT = timeit.default_timer()
CURR_ONE_HOUR_REQ_AT = timeit.default_timer()


def get_auth(payload):

    '''Get the token to start requesting data'''

    try:
        response = requests.post(
            'https://api.netatmo.com/oauth2/token', data=payload)

        response.raise_for_status()

        return response.json()

    except requests.exceptions.HTTPError as error:
        print(error.response.status_code, error.response.text)
    return


def refresh_auth(payload):

    '''Refresh the token'''

    try:
        response = requests.post(
            'https://api.netatmo.com/oauth2/token', data=payload)

        response.raise_for_status()

        return response.json()

    except requests.exceptions.HTTPError as error:
        print(error.response.status_code, error.response.text)
    return


def refresh_token(payload_auth, response_auth):

    '''Refresh a given token if its time is up'''

    global LAST_REFRESH_AUTH_AT, AUTH_EXP_IN

    if LAST_REFRESH_AUTH_AT is None:
        LAST_REFRESH_AUTH_AT = 0

    if 'expires_in' in response_auth:
        AUTH_EXP_IN = int(response_auth['expires_in'])

    elif 'expire_in' in response_auth:
        AUTH_EXP_IN = int(response_auth['expire_in'])

    else:
        raise RuntimeError(
            'Could not find when token expires: %s' % str(response_auth))

    LAST_REFRESH_AUTH_AT = timeit.default_timer()

    if LAST_REFRESH_AUTH_AT <= (0.80 * AUTH_EXP_IN):
        return

    print('\n' * 3)
    print('Refreshing Token...')

    payload_refresh = {
        'grant_type': 'refresh_token',
        'client_id': payload_auth['client_id'],
        'client_secret': payload_auth['client_secret'],
        'refresh_token': response_auth['refresh_token'],
        'scope': 'read_station'}

    response_refresh = refresh_auth(payload_refresh)

    response_auth['access_token'] = response_refresh['access_token']
    response_auth['refresh_token'] = response_refresh['refresh_token']
    LAST_REFRESH_AUTH_AT = 0
    print('\n' * 3)
    return


def check_reqs_limit(payload_auth, response_auth):

    '''To limit requests under 50 per 10 secs and 500 in an hour'''

    global TEN_SEC_REQ_CTR, ONE_HOUR_REQ_CTR
    global LAST_TEN_SEC_REQ_AT, CURR_TEN_SEC_REQ_AT
    global LAST_ONE_HOUR_REQ_AT, CURR_ONE_HOUR_REQ_AT

    if ((TEN_SEC_REQ_CTR > 40) and
        ((CURR_TEN_SEC_REQ_AT - LAST_TEN_SEC_REQ_AT) < 8)):

        TEN_SEC_REQ_CTR = 0

        print(
            'Too many reqs in 10 secs. Waiting for 15 secs. Time is:',
            time.asctime())

        time.sleep(15)
        LAST_TEN_SEC_REQ_AT = timeit.default_timer()

    if ((ONE_HOUR_REQ_CTR > 450) and
        ((CURR_ONE_HOUR_REQ_AT - LAST_ONE_HOUR_REQ_AT) < 3500)):

        ONE_HOUR_REQ_CTR = 0
        print(
            'Too many reqs in 1 hour. Waiting for 1.1 hours. Time is:',
            time.asctime())

        time.sleep(1.1 * 3600)
        LAST_ONE_HOUR_REQ_AT = timeit.default_timer()

    if (CURR_TEN_SEC_REQ_AT - LAST_TEN_SEC_REQ_AT) > 10:
        LAST_TEN_SEC_REQ_AT = timeit.default_timer()

    if (CURR_ONE_HOUR_REQ_AT - LAST_ONE_HOUR_REQ_AT) > 3600:
        LAST_ONE_HOUR_REQ_AT = timeit.default_timer()

    refresh_token(payload_auth, response_auth)

    TEN_SEC_REQ_CTR += 1
    CURR_TEN_SEC_REQ_AT = timeit.default_timer()

    ONE_HOUR_REQ_CTR += 1
    CURR_ONE_HOUR_REQ_AT = timeit.default_timer()
    return


def get_public_data(public_data_params, payload_auth, response_auth):

    '''Get station data within given bounds'''

    check_reqs_limit(payload_auth, response_auth)

    try:
        response = requests.post(
            'https://api.netatmo.com/api/getpublicdata',
            params=public_data_params)

        response.raise_for_status()

        meta_data = response.json()['body']

        return meta_data

    except requests.exceptions.HTTPError as error:
        print(error.response.status_code, error.response.text)
    return


def get_measured_data(measure_data_params, payload_auth, response_auth):

    '''Get station data for a given parameter'''

    check_reqs_limit(payload_auth, response_auth)

    try:
        response = requests.post(
            'https://api.netatmo.com/api/getmeasure',
            params=measure_data_params)

        response.raise_for_status()

        measure_data = response.json()['body']

        return measure_data

    except requests.exceptions.HTTPError as error:
        print(error.response.status_code, error.response.text)
    return


def assemble_stn_measure_data(measure_data, stn_id):

    '''Convert station data to pandas time series

    tested for rainfall data only
    '''

    time_stamps_list = []
    values_list = []

    for _dict in measure_data:
        _beg_time = _dict['beg_time']

        try:
            _step_time = _dict['step_time']
        except KeyError:
            _step_time = 0
            assert len(_dict['value']) == 1, (
                'Step_time is zero and values are more than 1!')

        _values = _dict['value']

        for i, _value in enumerate(_values):
            time_stamps_list.append(_beg_time + ((i + 1) * _step_time))
            values_list.append(_value[0])

    curr_idx = pd.np.array(time_stamps_list, dtype=pd.np.int64)
    curr_vals = pd.np.array(values_list, dtype=float)

    curr_stn_df = pd.DataFrame(
        index=curr_idx,
        data=curr_vals,
        columns=[stn_id],
        dtype=float)

    curr_stn_df.sort_index(inplace=True)
    curr_stn_df.index = pd.to_datetime(curr_stn_df.index, unit='s')
    return curr_stn_df
