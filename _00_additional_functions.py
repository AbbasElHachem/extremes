# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
GOAL: A combination of differrent functions used along the different scripts

Functions are called in different scripts, refer to the script to know where 
    it is called and what modifications may be needed
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# ============================================================================
import os

import fnmatch
import pyproj
import shapefile
import fiona
import osr

import numpy as np
import pandas as pd
import seaborn as sn
# import wradlib as wrl

import scipy.spatial as spatial
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# import scipy.stats
from matplotlib import path
from adjustText import adjust_text
from scipy.optimize import curve_fit
from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

# TODO: work more on documentation
#==============================================================================
#
#==============================================================================


def list_all_full_path(ext, file_dir):
    """
    Purpose: To return full path of files in all dirs of a given folder with a
    -------  given extension in ascending order.

    Keyword arguments:
    ------------------
        ext (string) = Extension of the files to list
            e.g. '.txt', '.tif'.
        file_dir (string) = Full path of the folder in which the files
            reside.
    """
    new_list = []
    patt = '*' + ext
    for root, _, files in os.walk(file_dir):
        for elm in files:
            if fnmatch.fnmatch(elm, patt):
                full_path = os.path.join(root, elm)
                new_list.append(full_path)
    return(sorted(new_list))

#==============================================================================
#
#==============================================================================


def convert_coords_fr_wgs84_to_utm32_(epgs_initial_str, epsg_final_str,
                                      first_coord, second_coord):
    """
    Purpose: Convert points from one reference system to a second
    --------
        In our case the function is used to transform WGS84 to UTM32
        (or vice versa), for transforming the DWD and Netatmo station
        coordinates to same reference system.

        Used for calculating the distance matrix between stations

    Keyword argument:
    -----------------
        epsg_initial_str: EPSG code as string for initial reference system
        epsg_final_str: EPSG code as string for final reference system
        first_coord: numpy array of X or Longitude coordinates
        second_coord: numpy array of Y or Latitude coordinates

    Returns:
    -------
        x, y: two numpy arrays containing the transformed coordinates in 
        the final coordinates system
    """
    initial_epsg = pyproj.Proj(epgs_initial_str)
    final_epsg = pyproj.Proj(epsg_final_str)
    x, y = pyproj.transform(initial_epsg, final_epsg,
                            first_coord, second_coord)
    return x, y

#==============================================================================
# Get colors by name
#==============================================================================


colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Sort colors by hue, saturation, value and name.
by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                for name, color in colors.items())
sorted_names_all = [name for hsv, name in by_hsv]

# for every 5min time period define a color, an ugly way of doing it, but easy


def pltcolor():
    cols_time_dict = {i: [] for i in np.arange(-60, 61, 5)}

    for i, l in enumerate(cols_time_dict.keys()):
        if l == -60:
            cols_time_dict[l].append(sorted_names_all[i * 12])
        elif l == -55:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -50:
            cols_time_dict[l].append(sorted_names_all[i * 8])
        elif l == -45:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -40:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -35:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -30:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -25:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -20:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -15:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -10:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -5:
            cols_time_dict[l].append(sorted_names_all[i * 6])

        elif l == 0:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 5:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 10:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 15:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 20:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 25:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 30:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 35:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 40:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 45:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 50:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 55:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 60:
            cols_time_dict[l].append(sorted_names_all[i * 6])
    return cols_time_dict

#==============================================================================
#
#==============================================================================


def get_stns_ids_coords_file(netatmo_coords_df_file):
    """
    Purpose: Get all Netatmo station Ids from coordinates dataframe
    -------
        Read the coordinates dataframe of the Netatmo data containing
        station Ids as index, convert each id to it'S original format
        (replace '_' with ':'), make it consistent with the original data

    Keyword arguments:
    ------------------
        netatmo:coords_df_file: path to the Netatmo coordinates dataframe
        make sure the station Ids are in the first column (index=0) and the
        data are seperated by a semi-column ';'.

    Returns:
    --------
        List of all station Ids in their original form as found in Stn_Mac
    """
    in_coords_df = pd.read_csv(netatmo_coords_df_file, index_col=0, sep=';')
    ids_lst = [stn_id.replace(':', '_') for stn_id in in_coords_df.index]
    ids_lst_unique = list(set(ids_lst))
    return ids_lst_unique
# @profile

#==============================================================================
#
#==============================================================================


def split_df_file_to_get_alls_stn_ids(all_df_files_list):
    """
    Purpose: Get Netatmo station Id from Netatmo data file path 
    -------
        Split the Netamo file path to get the station Id.
        Test to make sure it works

    Keyword arguments:
    ------------------
        all_df_files_list: a list of files (usually from list_all_full_path)
        each file refers to the path of a Netatmo data file

    Returns:
    --------
        List of all UNIQUE Ids of all stations found in the given Files 

    """
    # TODO: add example of file
    stns_ids = []
    for df_file in all_df_files_list:
        try:
            _, _, p3 = df_file.split('(')
            p3 = p3.replace(')).csv', '')
        except Exception as msg:
            print(msg)
        stns_ids.append(p3)
    stns_ids_unique_lst = list(set(stns_ids))
    return stns_ids_unique_lst

#==============================================================================
#
#==============================================================================


def split_one_df_file_to_get_stn_id(df_file):
    """
    Purpose: Get Station Id from Netatmo file path
    -------
    Keyword arguments:
    ------------------
        df_file: complete path the the Netatmo data file
    Return:
    -------
        Netatmo station ID
    """
    try:
        _, _, p3 = df_file.split('(')
        p3 = p3.replace(')).csv', '')
    except Exception as msg:
        print(msg)
    return p3

#==============================================================================
#
#==============================================================================


def ids_list_intersection(lst1, lst2):
    """
    Purpose: Intersect two lists
    -------
        Given two lists, find the intersection and save to a new list
    Keyword arguments:
    ------------------
        lst1: First list
        lst2: Second list
    Returns:
    --------
        New list containing the intersection of the two previous lists
    """
    # Use of hybrid method
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3

#==============================================================================
#
#==============================================================================


def remove_one_elt_lst(elt, lst):
    """ 
    Purpose: Remove a specific element from a list
    -------
        Note: removes all occurrences of the element
    Keyword arguments:
        elt: the element to remove (could be string or integer)
        lst: list containing all the different elements
    Returns:
        New list without the given element

    """
    new_lst = [el for el in lst if el != elt]
    return new_lst

#==============================================================================
#
#==============================================================================


def resampleDf(df, agg, closed='right', label='right',
               shift=False, leave_nan=True,
               label_shift=None,
               temp_shift=0,
               max_nan=0):
    """
    Purpose: Aggregate precipitation data

    Parameters:
    -----------
    Df: Pandas DataFrame Object
        Data set
    agg: string
        Aggregation 'M':Monthly 'D': Daily, 'H': Hourly, 'Min': Minutely
    closed: string
        'left' or 'right' defines the aggregation interval
    label: string
        'left' or 'right' defines the related timestamp
    shift: boolean, optional
        Shift the values by 6 hours according to the dwd daily station.
        Only valid for aggregation into daily aggregations
        True, data is aggregated from 06:00 - 06:00
        False, data is aggregated from 00:00 - 00:00
        Default is False

    temp_shift: shift the data based on timestamps (+- 0 to 5), default: 0

    label_shift: shift time label by certain values (used for timezones)

    leave_nan: boolean, optional
        True, if the nan values should remain in the aggregated data set.
        False, if the nan values should be treated as zero values for the
        aggregation. Default is True

    Remark:
    -------
        If the timestamp is at the end of the timeperiod:

        Input: daily        Output: daily+
            >> closed='right', label='right'

        Input: subdaily     Output: subdaily
            >> closed='right', label='right'

        Input: subdaily     Output: daily
            >> closed='right', label='left'

        Input: subdaily     Output: monthly
            >> closed='right', label='right'


        ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
        ! ! Always check, if aggregation is correct ! !
        ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !



    """

    if shift == True:
        df_copy = df.copy()
        if agg != 'D' or agg != '1440min':
            raise Exception('Shift can only be applied to daily aggregations')
        df = df.shift(-6, 'H')

    # To respect the nan values
    if leave_nan == True:
        # for max_nan == 0, the code runs faster if implemented as follows
        if max_nan == 0:
            print('Resampling')
            # Fill the nan values with values very great negative values and later
            # get the out again, if the sum is still negative
            df = df.fillna(-100000000000.)
            df_agg = df.resample(agg,
                                 closed=closed,
                                 label=label,
                                 base=temp_shift,
                                 loffset=label_shift).sum()
            # Replace negative values with nan values
            df_agg.values[df_agg.values[:] < 0.] = np.nan
        else:
            df_agg = df.resample(rule=agg,
                                 closed=closed,
                                 label=label,
                                 base=temp_shift,
                                 loffset=label_shift).sum()
            # find data with nan in original aggregation
            g_agg = df.groupby(pd.Grouper(freq=agg,
                                          closed=closed,
                                          label=label))
            n_nan_agg = g_agg.aggregate(lambda x: pd.isnull(x).sum())

            # set aggregated data to nan if more than max_nan values occur in the
            # data to be aggregated
            filter_nan = (n_nan_agg > max_nan)
            df_agg[filter_nan] = np.nan

    elif leave_nan == False:
        df_agg = df.resample(agg,
                             closed=closed,
                             label=label,
                             base=temp_shift,
                             loffset=label_shift).sum()
    if shift == True:
        df = df_copy
    return df_agg
#==============================================================================
#
#==============================================================================


def resample_Humidity_Df(data_frame,
                         temp_freq,
                         temp_shift=0,
                         label_shift=None,
                         df_sep_=None,
                         out_save_dir=None,
                         fillnan=False,
                         df_save_name=None,
                         method='mean'):
    """
    Purpose: Sample DF based on freq and time shift and label shift
    -------

    Keyword arguments: 
        data_frame: Netatmo humidity data frame to resample
        temp_freq: to which temporal frequency dataframe should be resampled
        temp_shift: shift the data based on timestamps (+- 0 to 5), default: 0
        label_shift: shift time label by certain values (used for timezones)
        df_sep_: seperater of the dataframe, used when saving
        out_save_dir: output location for saved dataframe
        fillnan: if True, fill nans with a certain value (usually 0 or -999)
        df_save_name: output name of resampled df
        method: which aggregation method to use when resampling (default: mean)
            select either the average or maximum or minimum value of aggregation
            period, could mean big difference for Humidity

    Returns:
    -------
        Resampled Netatmo Hmidity dataframe

    """
    # df_ = data_frame.copy()
    if method == 'mean':
        df_res = data_frame.resample(temp_freq,
                                     label='right',
                                     closed='right',
                                     loffset=label_shift,
                                     base=temp_shift).mean()
    if method == 'max':
        df_res = data_frame.resample(temp_freq,
                                     label='right',
                                     closed='right',
                                     loffset=label_shift,
                                     base=temp_shift).max()
    if method == 'min':
        df_res = data_frame.resample(temp_freq,
                                     label='right',
                                     closed='right',
                                     loffset=label_shift,
                                     base=temp_shift).min()
    if fillnan:
        df_res.fillna(value=0, inplace=True)
    if df_save_name is not None and out_save_dir is not None:
        df_res.to_csv(os.path.join(out_save_dir, df_save_name),
                      sep=df_sep_)
    return df_res


#==============================================================================
#
#==============================================================================

def select_df_within_period(df,  # original dataframe
                            start,  # startdate of period
                            end  # enddate of period
                            ):
    """ 
    Purpose: a function to select DF between two dates
    --------

    Keyword arguments:
    ------------------
        df: original dataframe
        start: start date of wanted period (should be same format as df index)
        end: end date of wanted period (should be same format as df index)
    Return:
    ------
        DF between giver start and end date (could return also empty Dataframe)
    """
    mask = (df.index > start) & (df.index <= end)
    df_period = df.loc[mask]
    return df_period

#==============================================================================
#
#==============================================================================


def select_df_within_period_year_basis(df,  # original dataframe
                                       start_day_of_year=91,  # start day april
                                       end_day_of_year=288  # end day 15oktober
                                       ):
    """ 
    Purpose: a function to select DF between two days, considering all years
    --------

    Keyword arguments:
    ------------------
        df: original dataframe
        start: start day of wanted period (should be same format as df index)
        end: end day of wanted period (should be same format as df index)
    Return:
    ------
        DF between giver start and end date (could return also empty Dataframe)
    """
    df = df.copy()

    mask = (df.index.dayofyear > start_day_of_year) & (
        df.index.dayofyear <= end_day_of_year)
    df_period = df.loc[mask]

    return df_period
#==============================================================================
#
#==============================================================================


def resample_intersect_2_dfs(df1,  # first dataframe to resample
                             df2,  # second dataframe to resample
                             temp_freq  # temporal frequency for resampling
                             ):
    """ 
    Purpose: Eesample and intersect two dataframes
    -------
        Also to handle Nans is Nans are found after the stations 
        have been intersected, nans are removed and reintersection is done

        Some workaround was done to deal with pandas Series and Dataframes
        as the Netatmo and DWD stations are different
    """
    df_resample1 = resampleDf(df=df1, agg=temp_freq)
    df_resample2 = resampleDf(df=df2, agg=temp_freq)
    print('\n+++ Done resampling each df, finding intersection+++')
    idx_common = df_resample1.index.intersection(df_resample2.index)

    if idx_common.shape[0] > 0:
        print('\n********\n common index is found, interescting dataframes')
        try:
            df_common1 = df_resample1.loc[idx_common, :]
            df_common2 = df_resample2.loc[idx_common, :]
        except Exception:

            df_common1 = df_resample1.loc[idx_common]
            df_common2 = df_resample2.loc[idx_common]
#             print(df_common2)
        df_common1 = pd.DataFrame(index=df_common1.index,
                                  data=df_common1.values)
        df_common2 = pd.DataFrame(index=df_common2.index,
                                  data=df_common2.values)
        print('\n********\n After resampling sum of NaN values is')
        print('first station has ', df_common1.isna().sum(), 'NaN values')
        print('second station has ', df_common2.isna().sum(), 'NaN values')

        try:
            if df_common1.isna().sum()[0] > 0:
                ix_df1_nans = df_common1.index[np.where(df_common1.isna())[0]]
                print('first station has Nans on ', len(ix_df1_nans), ' dates')
                df_common1.dropna(inplace=True)
        except IndexError:
            df_common1 = df_common1

        try:
            if df_common2.isna().sum()[0] > 0:
                ix_df2_nans = df_common2.index[np.where(df_common2.isna())[0]]
                print('second station has Nans on', len(ix_df2_nans), ' dates')
                df_common2.dropna(inplace=True)
        except IndexError:
            df_common2 = df_common2

        assert df_common1.isna().sum()[0] == 0, 'Nans in First station'
        assert df_common2.isna().sum()[0] == 0, 'Nans in second station'

        new_idx_common = df_common1.index.intersection(df_common2.index)
        if new_idx_common.shape[0] != idx_common.shape[0]:
            print('\n********\n Nan dropped fom Index, dfs resampled')

        try:
            df_common1 = df_common1.loc[new_idx_common, :]
            df_common2 = df_common2.loc[new_idx_common, :]
        except Exception:
            df_common1 = df_common1.loc[idx_common]
            df_common2 = df_common2.loc[idx_common]

        return df_common1, df_common2
    else:
        print('\n********\n no common index found, returning empty dfs')
        return pd.DataFrame(), pd.DataFrame()

#==============================================================================
#
#==============================================================================


def calculate_probab_ppt_below_thr(ppt_data, ppt_thr):
    """ calculate probability of values being below threshold """
    origin_count = ppt_data.shape[0]
    count_below_thr = ppt_data[ppt_data <= ppt_thr].shape[0]
    p0 = np.divide(count_below_thr, origin_count)
    return p0


#==============================================================================
#
#==============================================================================

# both correct
# def build_edf_fr_vals(ppt_data):
#     # Construct EDF, need to check if it works
#     """ construct empirical distribution function given data values """
#     data_sorted = np.sort(ppt_data, axis=0)[::-1]
#     x0 = np.squeeze(data_sorted)[::-1]
#     y0 = (np.arange(data_sorted.size) / len(data_sorted))
#     return x0, y0

#
def build_edf_fr_vals(data):
    """ construct empirical distribution function given data values """
    from statsmodels.distributions.empirical_distribution import ECDF
    data = data.ravel()
    cdf = ECDF(data)
    x0 = cdf.x[1:]
    y0 = cdf.y[1:]
    return x0, y0
#==============================================================================
#
#==============================================================================


def get_cdf_part_abv_thr(ppt_data, ppt_thr):
    """ select part of the CDF that is abv ppt thr """

    #p0 = calculate_probab_ppt_below_thr(ppt_data, ppt_thr)

    x0, y0 = build_edf_fr_vals(ppt_data)
    x_abv_thr = x0[x0 > ppt_thr]
    y_abv_thr = y0[np.where(x0 > ppt_thr)]

    #assert y_abv_thr[0] == p0, 'something is wrong with probability cal'

    return x_abv_thr, y_abv_thr
#==============================================================================
#
#==============================================================================


def get_dwd_stns_coords(coords_df_file, x_col_name, y_col_name, index_col_name,
                        sep_type):
    '''function used to return to coordinates dataframe of the DWD stations'''
    in_coords_df = pd.read_csv(coords_df_file, sep=sep_type,
                               dtype={index_col_name: str},
                               engine='c')
    in_coords_df.set_index(index_col_name, inplace=True)
    stn_ids = in_coords_df.index
    x_vals = in_coords_df[x_col_name].values.ravel()
    y_vals = in_coords_df[y_col_name].values.ravel()
    return in_coords_df, x_vals, y_vals, stn_ids
#==============================================================================
#
#==============================================================================


def get_nearest_dwd_station(first_stn_id, coords_df_file,
                            x_col_name, y_col_name, index_col_name,
                            sep_type, neighbor_to_chose):
    ''' Find for one station, the closest neibhouring station'''
    # read df coordinates and get station ids, and x, y values
    in_coords_df, x_vals, y_vals, stn_ids = get_dwd_stns_coords(
        coords_df_file, x_col_name, y_col_name, index_col_name, sep_type)
    # make a tupples from the coordinates
    coords_tuples = np.array([(x, y) for x, y in zip(x_vals, y_vals)])
    # create a tree from coordinates
    points_tree = spatial.cKDTree(coords_tuples)

    xstn = in_coords_df.loc[first_stn_id, x_col_name]
    ystn = in_coords_df.loc[first_stn_id, y_col_name]

    distances, indices = points_tree.query([xstn, ystn], k=10)
    coords_nearest_nbr = coords_tuples[indices[neighbor_to_chose]]
    stn_near = str(stn_ids[indices[neighbor_to_chose]])
    distance_near = distances[neighbor_to_chose]

    return coords_nearest_nbr, stn_near, distance_near

#==============================================================================
#
#==============================================================================


def get_netatmo_stns_coords(coords_df_file, x_col_name, y_col_name):
    """function used to return to coordinates dataframe of the DWD stations"""
    in_coords_df = pd.read_csv(coords_df_file, sep=';',
                               index_col=0, engine='c')
    stn_ids = list(set([stn_id.replace(':', '_')
                        for stn_id in in_coords_df.index]))

    lon_vals = in_coords_df[x_col_name].values.ravel()
    lat_vals = in_coords_df[y_col_name].values.ravel()

    return in_coords_df, lon_vals, lat_vals, stn_ids


#==============================================================================
#
#==============================================================================


def get_for_netatmo_nearest_dwd_station(first_stn_id, dwd_coords_df_file,
                                        netatmo_coords_df_file,
                                        x_col_name, y_col_name,
                                        lon_col_name, lat_col_name):
    """ Find for one netatmo station, the closest DWD neibhouring station"""

    wgs82 = "+init=EPSG:4326"
    utm32 = "+init=EPSG:32632"

    # read df coordinates and get station ids, and x, y values
    _, x_vals, y_vals, stn_ids = get_dwd_stns_coords(
        dwd_coords_df_file, x_col_name, y_col_name)

    in_netatmo_coords_df, _, _, _ = get_netatmo_stns_coords(
        netatmo_coords_df_file, lon_col_name, lat_col_name)

    first_stn_id = first_stn_id.replace('_', ':')
    lon_stn = in_netatmo_coords_df.loc[first_stn_id, lon_col_name]
    lat_stn = in_netatmo_coords_df.loc[first_stn_id, lat_col_name]

    x_stn, y_stn = convert_coords_fr_wgs84_to_utm32_(wgs82, utm32,
                                                     lon_stn, lat_stn)
    # make a tupples from the coordinates
    coords_tuples = np.array([(x, y) for x, y in zip(x_vals, y_vals)])

    # create a tree from coordinates
    points_tree = spatial.cKDTree(coords_tuples)

    distances, indices = points_tree.query([x_stn, y_stn], k=2)

    coords_nearest_nbr = coords_tuples[indices[1]]
    stn_near = str(stn_ids[indices[1]])
    distance_near = distances[1]

    return coords_nearest_nbr, stn_near, distance_near

#==============================================================================
#
#==============================================================================


def constrcut_contingency_table(stn1_id, stn2_id,
                                dataframe1, dataframe2, thr1, thr2):
    """ return percentage of values below or above a threshold
    return a matrix that is the contingency table"""
    df_1 = dataframe1.values.ravel().copy()
    df_2 = dataframe2.values.ravel().copy()

    assert df_1.shape[0] == df_2.shape[0], 'values should have same shape'

    df_combined = pd.DataFrame(index=dataframe1.index,
                               data=dataframe1.values,
                               columns=[stn1_id])

    # when two different stations (dwd-netatmo and not netatmo-netatmo)
    if stn1_id != stn2_id:
        df_combined.loc[:, stn2_id] = df_2
    else:
        stn2_id = stn2_id + '_'
        df_combined.loc[:, stn2_id] = df_2

    df_both_below_thr = np.nansum((df_combined[stn1_id].values <= thr1) & (
        df_combined[stn2_id].values <= thr2)) / df_1.shape[0]

    df_first_abv_second_below_thr = np.nansum((df_combined[stn1_id].values > thr1) & (
        df_combined[stn2_id].values <= thr2)) / df_1.shape[0]

    df_first_below_second_abv_thr = np.nansum((df_combined[stn1_id].values <= thr1) & (
        df_combined[stn2_id].values > thr2)) / df_1.shape[0]

    df_both_abv_thr = np.nansum((df_combined[stn1_id].values > thr1) & (
        df_combined[stn2_id].values > thr2)) / df_1.shape[0]

    return (100 * df_both_below_thr, 100 * df_first_abv_second_below_thr,
            100 * df_first_below_second_abv_thr, 100 * df_both_abv_thr)
#==============================================================================
#
#==============================================================================


def select_ppt_vals_abv_temp_thr(
    ppt_stn_id,  # netatmo ppt station id
    netatmo_ppt_stn1,  # df netatmo ppt station
    netatmo_temperature_df_file,  # path to df of all netatmo temp
    distance_matrix_netatmo_ppt_netatmo_temp,
    min_dist_thr_temp,  # distance threshold, selecting netatmo neigbours
    temp_thr,  # temp threshold, below assume ppt is snow
):
    """
    could be USED IN SCRIPT NUMBER  19, currently NOT USED
    """
    in_netatmo_temperature_stns_df = pd.read_csv(netatmo_temperature_df_file,
                                                 index_col=0, sep=';',
                                                 parse_dates=True,
                                                 infer_datetime_format=True,
                                                 engine='c')

    # read netatmo distance matrix ppt-temp
    in_df_distance_netatmo_netatmo = pd.read_csv(
        distance_matrix_netatmo_ppt_netatmo_temp, sep=';', index_col=0)
    # drop all duplicated stations in index
    in_df_distance_netatmo_netatmo = in_df_distance_netatmo_netatmo.loc[
        ~in_df_distance_netatmo_netatmo.index.duplicated(keep='first')]

    # Flag if consider Temperature Threshold
    # remove all days with Temp < temp_thr (1�C)

    # get distance netatmo ppt stn to all temp stns, get min
    distances_netatmo_to_stn1 = in_df_distance_netatmo_netatmo.loc[
        ppt_stn_id, :]
    sorted_distances_ppt_temp = distances_netatmo_to_stn1.sort_values(
        ascending=True)
    min_dist_ppt_temp = sorted_distances_ppt_temp.values[0]

    # check if temp station is close enough (5Km)
    if min_dist_ppt_temp <= min_dist_thr_temp:
        print('\n********\n Seperating distance is ',
              min_dist_ppt_temp, ' m')

        # if station is near, read df_temp
        netatmo_temp_stn = sorted_distances_ppt_temp.index[0]
        df_temp = in_netatmo_temperature_stns_df.loc[
            :, netatmo_temp_stn]
        df_temp.dropna(axis=0, inplace=True)

        print('\n********\n Temp Stn Id is', netatmo_temp_stn)
        # intersect netatmo ppt and netatmo temp data
        idx_cmn = netatmo_ppt_stn1.index.intersection(
            df_temp.index)
        df_netatmo_ppt_cmn = netatmo_ppt_stn1.loc[idx_cmn]
        df_temp_cmn = df_temp.loc[idx_cmn]

        # if enough data available, remove all temp vals < thr
        if (df_netatmo_ppt_cmn.values.shape[0] > 0 and
                df_temp_cmn.values.shape[0] > 0):
            df_ppt_abv_temp_thr = df_netatmo_ppt_cmn
            try:
                # find all dates where temp <= thr
                nan_vals_idx = df_temp_cmn[
                    df_temp_cmn <= temp_thr].index
                # -999
                df_ppt_abv_temp_thr.loc[nan_vals_idx] = np.nan

            except Exception as msg:
                print('error while selecting data abv thr', msg)

            else:
                print('station is near but no common data')

        else:
            print('\n********\n no nearby Temp station')

        if isinstance(df_ppt_abv_temp_thr, pd.Series):
            print('\n********\n Hours where temp<=1C ',
                  np.round(100 - 100 * df_ppt_abv_temp_thr.shape[0] /
                           df_netatmo_ppt_cmn.shape[0]), '%')
            netatmo_ppt_stn1 = df_ppt_abv_temp_thr

#==============================================================================
#
#==============================================================================


def select_convective_season(
    df,  # df to slice, index should be datetime
    month_lst  # list of month for convective season
):
    """
    return dataframe without the data corresponding to the winter season
    """
    df = df.copy()
    df_conv_season = df[~df.index.month.isin(month_lst)]

    return df_conv_season
#==============================================================================
#
#==============================================================================


def select_season(df,  # df to slice, index should be datetime
                  month_lst  # list of month for convective season
                  ):
    """
    return dataframe with the data corresponding to the season
    """
    df = df.copy()
    df_conv_season = df[df.index.month.isin(month_lst)]

    return df_conv_season

#==============================================================================
#
#==============================================================================


def compare_p1_dwd_p1_netatmo(
    val_dwd,  # p1 or ppt_mean  for df_DWD
    val_netatmo  # p1 or ppt_mean for df_Netatmo
):
    """
    use this function to compare two values from two different stations
    find if the values are equal +-10% , if smaller or if bigger
    return the result for being saved in a dataframe
    plot the result as scatter plots

    Return:
        marker, color 
    """
    _10_percent_dwd = 0.1 * val_dwd

    if val_dwd - _10_percent_dwd <= val_netatmo <= val_dwd + _10_percent_dwd:
        return 's', 'b'  # '='
    if val_netatmo < val_dwd - _10_percent_dwd:
        return '_', 'g'  # '-'
    if val_dwd + _10_percent_dwd < val_netatmo:
        return '+', 'r'  # '+'
    return

#==============================================================================
#
#==============================================================================


def look_agreement_df1_df2(stn_dwd_id,  # id of dwd station
                           stn_netatmo_id,  # if of netatmo station
                           df_dwd,  # df dwd intersected with netatmo
                           df_netatmo,  # df netatmo intersected with dwd
                           ppt_thr,  # ppt thr to check agreement
                           ):
    """
    read two dataframes, one for dwd and one for netatmo
    replace all values above threshold with 1 and below with 0

    calculate coorelation between the two
    look for agreement or disagreement
    try it with or without considering temperature threshold
    """
    df_dwd_copy = df_dwd.copy()
    df_netatmo_copy = df_netatmo.copy()

    # calculate pearson and spearman between original values
    abv_per_pear_corr = np.round(
        pears(df_dwd_copy.values,
              df_netatmo_copy.values)[0], 2)

    abv_per_spr_corr = np.round(spr(df_dwd_copy.values,
                                    df_netatmo_copy.values)[0], 2)

    # find nbr of events above thr
    events_dwd_abv_thr = df_dwd_copy[df_dwd_copy > ppt_thr].shape[0]
    events_netatmo_abv_thr = df_netatmo_copy[df_netatmo_copy >
                                             ppt_thr].shape[0]
    ratio_netatmo_dwd = np.round(
        (events_netatmo_abv_thr / events_dwd_abv_thr) * 100, 0)

    return (abv_per_pear_corr, abv_per_spr_corr,
            events_dwd_abv_thr, events_netatmo_abv_thr, ratio_netatmo_dwd)
#==============================================================================
#
#==============================================================================


def plot_subplot_fig(df_to_plot,    # dataframe to plot, coords-vals
                     col_var_to_plot,  # name of variable to plot, Marker value
                     col_color_of_marker,  # color of marker to use
                     plot_title,  # title of plot
                     lon_col_name='lon',  # name of longitude column
                     lat_col_name='lat',  # name of latitude column
                     shapefile_path=None,  # if Tur, plot Shapefile
                     table_to_plot=None,  # table with statistics, add to plot
                     out_save_name=None,  # Name of out figure
                     out_dir=None,  # out dir to use for saving
                     ):
    """
        Read the df_results containing for every netatmo station
        the coordinates (lon, lat) and the comparision between
        the p1 and the mean value, between netatmo and nearest dwd

        plot the reults on a map, either with or without temp_thr
        use the shapefile of BW
    """

    plt.ioff()

    gridsize = (4, 2)
    fig = plt.figure(figsize=(16, 12), dpi=150)
    fig.subplots_adjust(wspace=0.005, hspace=0.005)
    ax = plt.subplot2grid(gridsize, (0, 0), colspan=2, rowspan=3)
    ax2 = plt.subplot2grid(gridsize, (3, 0), colspan=2, rowspan=1)

    gs1 = gridspec.GridSpec(gridsize[0], gridsize[1])
    gs1.update(wspace=0.005, hspace=0.005)
    # set the spacing between axes.

    if shapefile_path is not None:
        # read and plot shapefile (BW or Germany) should be lon lat
        shp_de = shapefile.Reader(shapefile_path)
        for shape_ in shp_de.shapeRecords():
            lon = [i[0] for i in shape_.shape.points[:][::-1]]
            lat = [i[1] for i in shape_.shape.points[:][::-1]]
            ax.scatter(lon, lat, marker='.', c='lightgrey', alpha=0.25, s=2)

    # plot the stations in shapefile, look at the results of p1 comparasion
    for i in range(df_to_plot.shape[0]):
        ax.scatter(df_to_plot.loc[:, lon_col_name].values[i],
                   df_to_plot.loc[:, lat_col_name].values[i],
                   marker=df_to_plot.loc[:, col_var_to_plot].values[i],
                   c=df_to_plot.loc[:, col_color_of_marker].values[i],
                   alpha=1,
                   s=15,
                   label=str(col_var_to_plot))

    if table_to_plot is not None:
        table_to_plot = table_to_plot.T
        columns_to_plot = []
        for table_col in table_to_plot.columns:
            if col_var_to_plot in table_col:
                columns_to_plot.append(table_col)
        assert len(columns_to_plot) > 0, 'no columns found for table'
        table_col_to_plot = table_to_plot.loc[:, columns_to_plot]
        ax2.axis('tight')
        ax2.axis('off')
        table_plot = ax2.table(cellText=np.array([table_col_to_plot.values[0]]),
                               rowLabels=table_col_to_plot.index,
                               cellLoc='right',
                               colWidths=[0.1] * 3,
                               colColours=['b', 'g', 'r'],
                               colLabels=table_col_to_plot.columns,
                               clip_on=True,
                               in_layout=True,
                               snap=True,
                               loc='center')
        table_plot.auto_set_font_size(False)
        table_plot.set_fontsize(10)
        table_plot.scale(2, 2)

    ax.grid(alpha=0.5)
    ax.set_title(plot_title)

#     ax.set_aspect(1.0)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    fig.tight_layout()
    plt.savefig(os.path.join(out_dir, out_save_name),
                frameon=True, papertype='a4',
                bbox_inches='tight', pad_inches=.2)
    plt.close('all')

    return df_to_plot

#==============================================================================
#
#==============================================================================


def calc_plot_contingency_table_2_stations(
    first_stn_id,  # id of first stations
    second_stn_id,  # id of second station
    first_df,  # df of first station
    second_df,  # df of second station
    first_thr,  # threshhold for first df
    second_thr,  # threshhold for second df
    seperating_dist,  # distance between 2 stations
    temp_freq,  # temporal resolution of data
    df_append_resutl,  # df containing all stations as idx, append result
    out_plot_dir,  # out save dir for plot
    plot_figures,  # flag is True plot all plots
):
    """
    construct and save contingency table between two stations
    usually a DWD and a Netatmo Station used in script number _20_
    """
    # construct continegency table
    (df_both_below_thr,
     df_first_abv_second_below_thr,
     df_first_below_second_abv_thr,
     df_both_abv_thr) = constrcut_contingency_table(
        first_stn_id,
        second_stn_id,
        first_df,
        second_df,
        first_thr,
        second_thr)

    conf_arr = np.array([[df_both_below_thr,
                          df_first_abv_second_below_thr],
                         [df_first_below_second_abv_thr,
                          df_both_abv_thr]])
    #======================================================
    if conf_arr.shape[0] > 0:
        print('\n++++ Filling DF contingency table++++\n')
        # append to df, used later for analysis,
        df_append_resutl.loc[
            first_stn_id, 'DWD neighbour'
        ] = second_stn_id

        df_append_resutl.loc[
            first_stn_id, 'Both below %0.1fmm'
            % first_thr
        ] = np.round(df_both_below_thr, 2)

        df_append_resutl.loc[
            first_stn_id,
            'Netatmo below DWD'
        ] = np.round(df_first_below_second_abv_thr, 2)

        df_append_resutl.loc[
            first_stn_id,
            'Netatmo above DWD'
        ] = np.round(df_first_abv_second_below_thr, 2)

        df_append_resutl.loc[
            first_stn_id, 'Both above %0.1fmm'
            % first_thr
        ] = np.round(df_both_abv_thr, 2)

        if plot_figures:
            print('\n++++ plotting contingency table++++\n')
            plt.ioff()
            fig = plt.figure(figsize=(16, 12), dpi=100)
            ax = fig.add_subplot(111)
            ax2 = ax.twinx()
            ax.set_aspect(1)

            _ = sn.heatmap(conf_arr, annot=True,
                           vmin=0.0, cmap=plt.get_cmap('jet'),
                           vmax=100.0, fmt='.2f')

            ax.set_title(
                'Contingency table \n '
                'Rainfall Threshold %0.1fmm '
                '%s vs %s \n distance: %0.1f m, time freq: %s\n'
                % (first_thr,
                   first_stn_id, second_stn_id,
                   seperating_dist, temp_freq))

            ax.set_yticks([0.25, 1.25])
            ax.set_yticklabels(
                ['Netatmo below and DWD above',
                    'Netatmo below and DWD below'],
                rotation=90, color='g')
            ax.set_xticks([])

            ax2.set_yticks([0.25, 1.25])
            ax2.set_yticklabels(
                ['Netatmo above and DWD below',
                 'Netatmo above and DWD above'],
                rotation=-90, color='m')
            plt.savefig(
                os.path.join(
                    out_plot_dir,
                    'contingency_table_as_a_sequence'
                    '_%s_%s_%s_data_abv_%0.1fmm_.png'
                    % (first_stn_id, second_stn_id,
                        temp_freq,
                       first_thr)))

            plt.close(fig)
            plt.clf()

        return df_append_resutl

#==============================================================================
#
#==============================================================================


def save_how_many_abv_same_below(
    df_p1_mean,  # df with netatmo stns and result of comparing
    temp_freq,  # temporal freq of df
    ppt_thr,  # used for calculating p1 or p5
    out_dir,  # out save dir for plots
    year_vals  # if all years or year by year
):
    """
    use this funcion to find how many stations have similar values for p1 and
    for the average, how many are netatmo are below dwd and how many 
    netatmo are above dwd

    return the result in a dataframe
    """
    # drop all netatmo stns with no data
    df_p1_mean.dropna(axis=0, how='any', inplace=True)

    df_out = pd.DataFrame()

    df_out.loc['count_vals_same_p%d' % ppt_thr,
               'Result'] = df_p1_mean[
                   df_p1_mean['p%d' % ppt_thr] == 's'].shape[0]
    df_out.loc['count_vals_less_p%d' % ppt_thr,
               'Result'] = df_p1_mean[
                   df_p1_mean['p%d' % ppt_thr] == '_'].shape[0]
    df_out.loc['count_vals_abv_p%d' % ppt_thr,
               'Result'] = df_p1_mean[
                   df_p1_mean['p%d' % ppt_thr] == '+'].shape[0]

    df_out.loc['count_vals_same_mean_ppt',
               'Result'] = df_p1_mean[
                   df_p1_mean['mean_ppt'] == 's'].shape[0]
    df_out.loc['count_vals_less_mean_ppt',
               'Result'] = df_p1_mean[
                   df_p1_mean['mean_ppt'] == '_'].shape[0]
    df_out.loc['count_vals_abv_mean_ppt',
               'Result'] = df_p1_mean[
                   df_p1_mean['mean_ppt'] == '+'].shape[0]

    df_out.to_csv(os.path.join(
        out_dir,
        'year_%s_df_similarities_%s_%dmm_.csv'
        % (year_vals, temp_freq, ppt_thr)))
    return df_out
#==============================================================================
#
#==============================================================================


def plt_on_map_agreements(
    df_correlations,  # df with netatmo stns, result of correlations
    col_to_plot,  # which column to plot , str_label_of_col
    shp_de_file,  # shapefile of BW
    temp_freq,  # temp freq of df
    out_dir,  # out save dir for plots
    year_vals,  # if all years or year by year
    val_thr_percent  # consider all values above it
):
    """
    Read the df_results containing for every netatmo station
    the coordinates (lon, lat) and the comparision between 
    the pearson and spearman correlations,
    between netatmo and nearest dwd station
    plot the reults on a map, either with or with temp_thr
    """
    if 'Bool' in col_to_plot:
        percent_add = '_above_%dpercent' % val_thr_percent
    else:
        percent_add = '_'
    print('plotting comparing %s' % col_to_plot)
    # select values only for column and dron nans
    df_correlations = df_correlations.loc[:, [
        'lon', 'lat', col_to_plot]].dropna()

    plt.ioff()
    fig = plt.figure(figsize=(15, 15), dpi=150)

    ax = fig.add_subplot(111)

    shp_de = shapefile.Reader(shp_de_file)
    # read and plot shapefile (BW or Germany) should be lon lat
    for shape_ in shp_de.shapeRecords():
        lon = [i[0] for i in shape_.shape.points[:][::-1]]
        lat = [i[1] for i in shape_.shape.points[:][::-1]]
        ax.scatter(lon, lat, marker='.', c='lightgrey',
                   alpha=0.25, s=2)

    # plot the stations in shapefile, look at the results of agreements
    texts = []
    for i in range(df_correlations.shape[0]):
        ax.scatter(df_correlations.lon.values[i],
                   df_correlations.lat.values[i],
                   alpha=1,
                   c='b',
                   s=15,
                   label=df_correlations[col_to_plot].values[i])
        texts.append(ax.text(df_correlations.lon.values[i],
                             df_correlations.lat.values[i],
                             str(df_correlations[col_to_plot].values[i]),
                             color='k'))

    ax.set_title('%s Data %s'
                 ' Netatmo and DWD %s data year %s'
                 % (col_to_plot,  percent_add,
                     temp_freq, year_vals))
    ax.grid(alpha=0.5)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect(1.0)
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='->', color='red', lw=0.25))

    plt.savefig(
        os.path.join(
            out_dir,
            'year_%s_%s_%s_netatmo_ppt_dwd_station_%s.png'
            % (year_vals, temp_freq, col_to_plot, percent_add)),
        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)
    plt.close()
    return df_correlations

#==============================================================================
#
#==============================================================================


def plt_correlation_with_distance(
    df_correlations,  # df with netatmo stns, result of correlations
    dist_col_to_plot,  # name of column with distance values
    corr_col_to_plot,  # which column to plot , str_label_of_col
    temp_freq,  # temp freq of df
    out_dir,  # out save dir for plots
    year_vals,  # if all years or year by year
    val_thr_percent,  # consider all values above it
    neighbor_nbr  # which neighbor was chosen
):
    """
    Read the df_results containing for every netatmo station
    the coordinates (lon, lat) and the comparision between 
    the pearson and spearman correlations,
    between netatmo and nearest dwd station
    plot the reults, correlation vs distance
    """
    if 'Bool' in corr_col_to_plot:
        percent_add = '_above_%dpercent' % val_thr_percent
    else:
        percent_add = '_'
    print('plotting comparing %s' % corr_col_to_plot)
    # select values only for column and dron nans
    df_correlations = df_correlations.loc[:, [
        'lon', 'lat', dist_col_to_plot, corr_col_to_plot]].dropna()

    plt.ioff()
    fig = plt.figure(figsize=(16, 12), dpi=150)

    ax = fig.add_subplot(111)

    # plot the stations in shapefile, look at the results of agreements
    colors_arr = (df_correlations.loc[:, corr_col_to_plot].values /
                  max(df_correlations.loc[:, corr_col_to_plot].values))
    colors_arr[colors_arr < 0] = 0
    ax.scatter(df_correlations.loc[:, dist_col_to_plot].values,
               df_correlations.loc[:, corr_col_to_plot].values,
               alpha=.8,
               c='b',  # colors_arr,
               s=15,
               marker='d',
               # cmap=plt.get_cmap('viridis'),
               label='Number of Stations %d' %
               df_correlations.loc[:, dist_col_to_plot].values.shape[0])

    ax.set_title('%s Data %s Neighbor number %d'
                 ' Netatmo and Netamto %s data year %s'
                 % (corr_col_to_plot,  percent_add, neighbor_nbr,
                     temp_freq, year_vals))
    ax.grid(alpha=0.25)
    ax.set_ylim([-0.1, 1.1])
    ax.set_xlim([0, max(df_correlations.loc[:, dist_col_to_plot].values)])
    ax.set_xlabel('Distance to neighbour (m)')
    ax.legend(loc=0)
    ax.set_ylabel('Indicator Spearman Correlation')

    plt.savefig(
        os.path.join(
            out_dir,
            'indic_corr_with_dist_%s_%s_%s_netatmo_ppt_netatmo_above_%s_neighbor_nbr_%d_.png'
            % (year_vals, temp_freq, corr_col_to_plot, percent_add, neighbor_nbr)),
        frameon=True, papertype='a4',
        bbox_inches='tight', pad_inches=.2)
    plt.close()
    return df_correlations

#==============================================================================
#
#==============================================================================


def func(x, a, b, c, d):
    """ 3degree polynomial function used for fitting and as filter"""
    return a * x**3 + b * x**2 + c * x + d


def fit_curve_get_vals_below_curve(x, y, func, stns, shift_per=0.15):
    """ fit function to data and shifted 10% downwards
    return all data above and below function with 
    corresponding distance and correlation values and station ids
    """

    # bounds=[[a1,b1],[a2,b2]]

    popt, _ = curve_fit(
        func, x, y)

    print('fitted parameters are ', popt[0], popt[1], popt[2])  # , popt[3])
    y_fitted = func(x, *popt)

    lower_bound = y_fitted.max() * shift_per
    y_fitted_shifted = y_fitted - lower_bound

    xvals_below_curve = x[np.where(y <= y_fitted_shifted)]
    yvals_below_curve = y[np.where(y <= y_fitted_shifted)]

    stns_below_curve = stns[np.where(y <= y_fitted_shifted)]

    xvals_above_curve = x[np.where(y > y_fitted_shifted)]
    yvals_above_curve = y[np.where(y > y_fitted_shifted)]

    stns_above_curve = stns[np.where(y > y_fitted_shifted)]

    return (y_fitted_shifted, xvals_below_curve, yvals_below_curve,
            xvals_above_curve, yvals_above_curve, stns_below_curve,
            stns_above_curve)


#==============================================================================
#
#==============================================================================
def gen_path_df_file(main_path, start_path_acc, time_freq,
                     data_source, percent, neighbor,
                     use_filtered_data=False,
                     filter_percent=None):
    """ use this funtion to get the path to the dfs for different neighbors"""

    if use_filtered_data is False:
        return main_path / (
            r'%sfreq_%s_%s_netatmo_upper_%s_percent_data_considered_neighbor_%d_.csv'
            % (start_path_acc, time_freq, data_source, percent, neighbor))

    else:
        return main_path / (
            r'%sfreq_%s_%s_netatmo_upper_%s_percent_data_considered_neighbor_%d_filtered_%s.csv'
            % (start_path_acc, time_freq, data_source, percent, neighbor, filter_percent))

#==============================================================================
#
#==============================================================================


def read_filter_df_corr_return_stns_x_y_vals(df_file, thr_low=0, thr_high=1,
                                             x_col_name='Distance to neighbor',
                                             y_col_name='Bool_Spearman_Correlation'):
    """ read df with correlation values, select all between 0 and 1 and 
        return station_ids, distance_values, correlation_values, and df
    """
    in_df = pd.read_csv(df_file, index_col=0, sep=';').dropna(how='any')
    # ycol = in_df.columns.intersection([y_col_name])
    # xcol = in_df.columns.intersection([x_col_name])
    in_df = in_df[(thr_low <= in_df.loc[:, y_col_name]) &
                  (in_df.loc[:, y_col_name] < thr_high)]
    stn_ids = in_df.index
    x_vals = in_df.loc[:, x_col_name].values.ravel()
    y_vals = in_df.loc[:, y_col_name].values.ravel()
    print('Number of stations is', len(stn_ids))
    return stn_ids, x_vals, y_vals, in_df
#==============================================================================
#
#==============================================================================


def remove_all_low_corr_short_dist(x, y, xthr, ythr, stns, df,
                                   x_col_name='Distance to neighbor',
                                   y_col_name='Bool_Spearman_Correlation'):
    """
    for all distances below x_thr remove the correlation values
    these are outliers that will affect the fitting of the curve
    get stations with high correlations and drop nans and reselect data
    return x_good, y_good, stns_good, df_good
    """
    df = df.copy()
    for i, dist in enumerate(x):
        if dist <= xthr:
            if y[i] <= ythr:
                y[i] = np.nan
                x[i] = np.nan
    # remove nans and get all stations and data
    s_gd_corr = stns[np.where(y > 0) and np.where(x > 0)]
    df_gd_corr = df.loc[s_gd_corr, :].dropna(how='any')
    x_gd_corr = df_gd_corr.loc[s_gd_corr, x_col_name].values.ravel()
    y_gd_corr = df_gd_corr.loc[s_gd_corr, y_col_name].values.ravel()

    # assert np.isnan(x_gd_corr).sum() == 0, 'nans in data'
    return x_gd_corr, y_gd_corr, s_gd_corr, df_gd_corr
#==============================================================================
#
#==============================================================================


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

#==============================================================================
#
#==============================================================================


def find_nearest(array, value):
    ''' given a value, find nearest one to it in original data array'''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
#==============================================================================
#
#==============================================================================


def get_radar_intense_events(radar_files_loc, intense_events_df_index_lst, temp_freq):
    ''' Given a list of radolan files and rainfall intense events df
        find the corresponding radolan files, save to a new df,
        index=Event Time; columns=file_path
    '''
    radolan_events_to_keep, event_dates = [], []
    for i, file in enumerate(radar_files_loc):
        #     file = file + '.gz'

        event_date_raw = file.split('\\')[-1].split('-')[2]
        event_date = ('20' + event_date_raw[:2] + '-' + event_date_raw[2:4] +
                      '-' + event_date_raw[4:6] + ' ' + event_date_raw[6:8] +
                      ':' + event_date_raw[8:] + ':00')

        event_date = pd.DatetimeIndex([event_date])

        event_end_date = event_date + pd.Timedelta(minutes=10)

        if event_end_date[0] in intense_events_df_index_lst:
            print('Event: ', i, '/', len(radar_files_loc),  event_date)
            radolan_events_to_keep.append(file)
            event_dates.append(event_end_date[0])
    df_radar_events_to_keep = pd.DataFrame(index=event_dates,
                                           data=radolan_events_to_keep)
    # TODO: change path
    df_radar_events_to_keep.to_csv(r'X:\exchange\ElHachem'
                                   r'\%s_intense_events_radolan_files.csv'
                                   % temp_freq,
                                   sep=';')

    return df_radar_events_to_keep

#==============================================================================
#
#==============================================================================


# def mask_radolan_based_on_shp(radolan_file, shp_file):
#     '''
#         read a radolan file, transform from Polar to WGS84
#         get longitude and latitude, create mask based on
#         shapefile location, save the mask using numpy dump
#
#     '''
#     rwdata, rwattrs = wrl.io.read_radolan_composite(radolan_file)
#     # mask data
#     sec = rwattrs['secondary']
#     rwdata.flat[sec] = -9999
#     rwdata = np.ma.masked_equal(rwdata, -9999)
#
#     # create radolan projection object
#     proj_stereo = wrl.georef.create_osr("dwd-radolan")
#
#     # create wgs84 projection object
#     proj_wgs = osr.SpatialReference()
#     proj_wgs.ImportFromEPSG(4326)
#
#     # get radolan grid
#     radolan_grid_xy = wrl.georef.get_radolan_grid(900, 900)
# #         # convert to lonlat
#     radolan_grid_ll = wrl.georef.reproject(radolan_grid_xy,
#                                            projection_source=proj_stereo,
#                                            projection_target=proj_wgs)
#     lon1 = radolan_grid_ll[:, :, 0]
#     lat1 = radolan_grid_ll[:, :, 1]
#     # read shapefile
#     ishape = fiona.open(shp_file)
#     first = ishape.next()
#
#     mask = np.ones_like(lon1, dtype=np.bool)
#
#     for _, i_poly in enumerate(first['geometry']['coordinates']):
#
#         p = path.Path(np.array(i_poly)[0, :, :])
#         grid_mask = p.contains_points(
#             np.vstack((lon1.flatten(),
#                        lat1.flatten())).T).reshape(900, 900)
#     mask[grid_mask] = 0
#     mask.dump(r'C:\Users\hachem\Desktop\radar\masked_array_bw')
#
#     return mask
