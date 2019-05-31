import os

import fnmatch
import pyproj

import matplotlib.colors as mcolors

import numpy as np
import pandas as pd
import scipy.spatial as spatial
#==============================================================================
#
#==============================================================================


def list_all_full_path(ext, file_dir):
    """
    Purpose: To return full path of files in all dirs of a 
            given folder with a
            given extension in ascending order.
    Description of the arguments:
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
    '''fct to convert points from wgs 84 to utm32'''
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
    ''' get all ids from coords_file'''
    in_coords_df = pd.read_csv(netatmo_coords_df_file, index_col=0, sep=';')
    ids_lst = [stn_id.replace(':', '_') for stn_id in in_coords_df.index]
    ids_lst_unique = list(set(ids_lst))
    return ids_lst_unique
# @profile

#==============================================================================
#
#==============================================================================


def split_df_file_to_get_alls_stn_ids(all_df_files_list):
    ''' get stn_id from netatmo file'''
    stns_ids = []
    for df_file in all_df_files_list:
        try:
            _, _, p3 = df_file.split('(')
            p3 = p3.replace(')).csv', '')
        except Exception as msg:
            print(msg)
        stns_ids.append(p3)
    stns_ids_unique = list(set(stns_ids))
    return stns_ids_unique

#==============================================================================
#
#==============================================================================


def split_one_df_file_to_get_stn_id(df_file):
    ''' get stn_id from netatmo file'''
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
    ''' intersect two lists'''
    # Use of hybrid method
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3

#==============================================================================
#
#==============================================================================


def remove_one_elt_lst(elt, lst):
    ''' remove a specific element from a list'''
    new_lst = []
    for el in lst:
        if el != elt:
            new_lst.append(el)
    return new_lst

#==============================================================================
#
#==============================================================================


def resampleDf(data_frame,
               temp_freq,
               temp_shift=0,
               label_shift=None,
               df_sep_=None,
               out_save_dir=None,
               fillnan=False,
               df_save_name=None):
    ''' sample DF based on freq and time shift and label shift '''

    # df_ = data_frame.copy()
    df_res = data_frame.resample(temp_freq,
                                 label='right',
                                 closed='right',
                                 loffset=label_shift,
                                 base=temp_shift).apply(lambda x: x.values.sum())

    if fillnan:
        df_res.fillna(value=0, inplace=True)
    if df_save_name is not None and out_save_dir is not None:
        df_res.to_csv(os.path.join(out_save_dir, df_save_name),
                      sep=df_sep_)
    return df_res

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
    ''' sample DF based on freq and time shift and label shift '''

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

def select_df_within_period(df, start, end):
    ''' a function to select df between two dates'''
    mask = (df.index > start) & (df.index <= end)
    df_period = df.loc[mask]
    return df_period

#==============================================================================
#
#==============================================================================


def resample_intersect_2_dfs(df1, df2, temp_freq):
    ''' a fct to resample and intersect two dataframes'''
    df_resample1 = resampleDf(data_frame=df1, temp_freq=temp_freq)
    df_resample2 = resampleDf(data_frame=df2, temp_freq=temp_freq)

    idx_common = df_resample1.index.intersection(df_resample2.index)

    try:
        df_common1 = df_resample1.loc[idx_common, :]
        df_common2 = df_resample2.loc[idx_common, :]
    except Exception:
        df_common1 = df_resample1.loc[idx_common]
        df_common2 = df_resample2.loc[idx_common]
    return df_common1, df_common2

#==============================================================================
#
#==============================================================================


def calculate_probab_ppt_below_thr(ppt_data, ppt_thr):
    ''' calculate probability of values being below threshold '''
    origin_count = ppt_data.shape[0]
    count_below_thr = ppt_data[ppt_data <= ppt_thr].shape[0]
    p0 = np.divide(count_below_thr, origin_count)
    return p0


#==============================================================================
#
#==============================================================================


def build_edf_fr_vals(ppt_data):
    # Construct EDF
    ''' construct empirical distribution function given data values '''
    data_sorted = np.sort(ppt_data, axis=0)[::-1]
    x0 = np.squeeze(data_sorted)[::-1]
    y0 = (np.arange(data_sorted.size) / len(data_sorted))
    return x0, y0


#==============================================================================
#
#==============================================================================


def get_cdf_part_abv_thr(ppt_data, ppt_thr):
    ''' select part of the CDF that is abv ppt thr '''

    p0 = calculate_probab_ppt_below_thr(ppt_data, ppt_thr)

    x0, y0 = build_edf_fr_vals(ppt_data)
    x_abv_thr = x0[x0 > ppt_thr]
    y_abv_thr = y0[np.where(x0 > ppt_thr)]
    assert y_abv_thr[0] == p0, 'something is wrong with probability cal'

    return x_abv_thr, y_abv_thr
#==============================================================================
#
#==============================================================================


def get_dwd_stns_coords(coords_df_file, x_col_name, y_col_name):
    '''function used to return to coordinates dataframe of the DWD stations'''
    in_coords_df = pd.read_csv(coords_df_file, sep=';',
                               index_col=3, engine='c')
    stn_ids = in_coords_df.index
    x_vals = in_coords_df[x_col_name].values.ravel()
    y_vals = in_coords_df[y_col_name].values.ravel()
    return in_coords_df, x_vals, y_vals, stn_ids

#==============================================================================
#
#==============================================================================


def get_netatmo_stns_coords(coords_df_file, x_col_name, y_col_name):
    '''function used to return to coordinates dataframe of the DWD stations'''
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
    ''' Find for one netatmo station, the closest DWD neibhouring station'''

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
    ''' return percentage of values below or above a threshold'''
    df_1 = dataframe1.values.ravel()
    df_2 = dataframe2.values.ravel()

    assert df_1.shape[0] == df_2.shape[0], 'values should have same shape'

    df_combined = dataframe1.copy()
    df_combined.loc[:, stn2_id] = df_2

    print(df_combined)

    df_both_below_thr = ((df_combined[stn1_id].values <= thr1) & (
        df_combined[stn2_id].values <= thr2)).sum() / df_1.shape[0]

    df_first_abv_second_below_thr = ((df_combined[stn1_id].values > thr1) & (
        df_combined[stn2_id].values <= thr2)).sum() / df_1.shape[0]

    df_first_below_second_abv_thr = ((df_combined[stn1_id].values <= thr1) & (
        df_combined[stn2_id].values > thr2)).sum() / df_1.shape[0]

    df_both_abv_thr = ((df_combined[stn1_id].values > thr1) & (
        df_combined[stn2_id].values > thr2)).sum() / df_1.shape[0]

    return (100 * df_both_below_thr, 100 * df_first_abv_second_below_thr,
            100 * df_first_below_second_abv_thr, 100 * df_both_abv_thr)
#==============================================================================
#
#==============================================================================
